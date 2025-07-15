#!/usr/bin/env python3
"""
Coherence Analyzer

A module for aggregating evidence and checking biological consistency across layers.
Performs tissue-specific biomarker relevance checks using GTEx API v2.
Supports decimal versions of GENCODE IDs (e.g., ENSG00000141510.15).

Example usage:
    analyzer = CoherenceAnalyzer()
    result = await analyzer.analyze_coherence(
        "processed_input.json",
        "examples/ex_clinical_data.json"
    )
"""

import os
import json
import logging
import asyncio
import hashlib
import time
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import statistics
import numpy as np

# Configure logging
from .logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)

# Async HTTP client
try:
    import aiohttp # type: ignore
    AIOHTTP_AVAILABLE = True
except ImportError:
    logger.warning("aiohttp not available. Install with: pip install aiohttp")
    AIOHTTP_AVAILABLE = False


@dataclass
class CoherenceResult:
    """Result of coherence analysis"""
    cross_scale_consistency: float
    tissue_specificity_scores: Dict[str, float] = field(default_factory=dict)
    evidence_weights: Dict[str, float] = field(default_factory=dict)
    gtex_validation_results: Dict[str, Any] = field(default_factory=dict)
    # warnings: List[str] = field(default_factory=list)
    processing_timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class CoherenceAnalyzer:
    """
    Aggregates evidence and checks biological consistency across layers
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the coherence analyzer
        
        Args:
            config: Configuration dictionary
        """
        if not AIOHTTP_AVAILABLE:
            raise ImportError("aiohttp is required for GTEx API calls")
        
        if config is None:
            logger.warning("No configuration provided, using default configuration")
            config = self._get_default_config()
        
        self.config = config
        self.cache_dir = Path("cache/gtex")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # API configurations
        self.gtex_api_base = "https://gtexportal.org/api/v2"
        self.session = None
        
        # Tissue mapping for disease contexts (using GTEx v2 tissue IDs)
        self.disease_tissue_mapping = {
            'breast_cancer': ['Breast_Mammary_Tissue'],
            'lung_cancer': ['Lung'],
            'liver_cancer': ['Liver'],
            'prostate_cancer': ['Prostate'],
            'colon_cancer': ['Colon_Sigmoid', 'Colon_Transverse'],
            'brain_cancer': ['Brain_Cortex', 'Brain_Frontal_Cortex_BA9'],
            'kidney_cancer': ['Kidney_Cortex'],
            'pancreatic_cancer': ['Pancreas'],
            'ovarian_cancer': ['Ovary'],
            'skin_cancer': ['Skin_Sun_Exposed_Lower_leg', 'Skin_Not_Sun_Exposed_Suprapubic'],
            'blood_cancer': ['Whole_Blood'],
            'stomach_cancer': ['Stomach']
        }
        
        logger.info("CoherenceAnalyzer initialized")
    
    def _get_default_config(self) -> Dict:
        """Get default configuration"""
        logger.warning("Using default configuration values")
        return {
            'evidence_weights': {
                'variant_consistency': 0.4,
                'tissue_specificity': 0.6
            },
            'gtex_cache_ttl': 604800,  # 7 days in seconds
            'gtex_timeout': 30,  # API timeout in seconds
            'gtex_max_retries': 1,
            'expression_threshold': 1.0,  # Minimum TPM for meaningful expression
            'tissue_specificity_threshold': 0.5,
            'consistency_threshold': 0.6
        }
     
    async def analyze_coherence(self, 
                              processed_input_file: str,
                              clinical_data_file: str) -> CoherenceResult:
        """
        Main analysis method for coherence across biological scales
        
        Args:
            processed_input_file: Path to processed input JSON
            clinical_data_file: Path to clinical data JSON
            
        Returns:
            CoherenceResult containing cross-scale consistency score
        """
        logger.info("Starting coherence analysis across biological scales")
        
        # Load input data
        processed_input = self._load_json(processed_input_file)
        clinical_data = self._load_json(clinical_data_file)
        
        # Initialize result object
        result = CoherenceResult(cross_scale_consistency=0.0)
        
        # Perform API-dependent analysis
        async with aiohttp.ClientSession(
            timeout=aiohttp.ClientTimeout(total=self.config['gtex_timeout'])
        ) as session:
            self.session = session

            # Extract genes and GENCODE IDs from processed input
            genes = self._extract_genes(processed_input)
            gene_to_gencode = self._extract_gencode_ids(processed_input)
            tissues = self._extract_target_tissues(clinical_data)
            
            # Get GTEx expression data using single gene-single tissue queries
            gtex_data = await self._get_gtex_expression_data_single_queries(genes, gene_to_gencode, tissues)
            result.gtex_validation_results = gtex_data
            
            # Calculate coherence across scales
            consistency_score = self._calculate_cross_scale_consistency(
                processed_input, clinical_data, gtex_data
            )
            
            result.cross_scale_consistency = consistency_score
            result.evidence_weights = self.config['evidence_weights']
            
            # Calculate tissue-specific scores
            result.tissue_specificity_scores = self._calculate_tissue_specificity_scores(
                clinical_data, gtex_data
            )
            
            # Add warnings if consistency is low
            # if consistency_score < self.config['consistency_threshold']:
            #     result.warnings.append(f"Low cross-scale consistency: {consistency_score:.3f}")
            
            logger.info(f"Coherence analysis complete. Consistency score: {consistency_score:.3f}")
        
        return result
    
    def _load_json(self, file_path: str) -> Dict:
        """Load JSON file with error handling"""
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            logger.error(f"File not found: {file_path}")
            return {}
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in {file_path}: {e}")
            return {}
    
    def _extract_genes(self, processed_input: Dict) -> List[str]:
        """Extract unique genes from processed input"""
        genes = set()
        
        if 'missense_variants' in processed_input:
            for variant in processed_input['missense_variants']:
                if 'gene' in variant and variant['gene']:
                    genes.add(variant['gene'])
        
        gene_list = list(genes)
        logger.info(f"Extracted {len(gene_list)} unique genes: {gene_list}")
        return gene_list
    
    def _extract_gencode_ids(self, processed_input: Dict) -> Dict[str, str]:
        """
        Extract GENCODE gene IDs directly from processed input.
        Supports decimal versions and implements version reduction logic.
        """
        gene_to_gencode = {}
        
        if 'missense_variants' in processed_input:
            for variant in processed_input['missense_variants']:
                gene = variant.get('gene')
                gencode_id = variant.get('gencode_id')
                
                if gene and gencode_id:
                    # Store the original GENCODE ID
                    gene_to_gencode[gene] = gencode_id
                    logger.debug(f"Extracted GENCODE ID for {gene}: {gencode_id}")

        logger.info(f"Extracted GENCODE IDs for {len(gene_to_gencode)} genes: {gene_to_gencode}")
        return gene_to_gencode

    def _generate_version_variants(self, gencode_id: str, max_attempts: int = 10) -> List[str]:
        """
        Generate version variants of a GENCODE ID by reducing version number.
        
        Args:
            gencode_id: Original GENCODE ID (e.g., "ENSG00000141510.19")
            max_attempts: Maximum number of version reductions to try
            
        Returns:
            List of GENCODE IDs to try, starting with original
        """
        if not gencode_id or '.' not in gencode_id:
            return [gencode_id]
        
        base_id, version_str = gencode_id.split('.', 1)
        try:
            version = int(version_str)
        except ValueError:
            logger.warning(f"Invalid version number in GENCODE ID: {gencode_id}")
            return [gencode_id]
        
        variants = [gencode_id]  # Start with original
        
        # Generate version variants by reducing version number
        for i in range(1, max_attempts + 1):
            new_version = version - i
            if new_version >= 0:  # Don't go below 0
                variants.append(f"{base_id}.{new_version}")
        
        logger.debug(f"Generated {len(variants)} version variants for {gencode_id}: {variants}")
        return variants

    def _extract_target_tissues(self, clinical_data: Dict) -> List[str]:
        """Extract relevant tissues based on clinical diagnosis"""
        diagnosis = clinical_data.get('diagnosis', '').lower()
        
        # Map diagnosis to GTEx v2 tissues
        tissues = self.disease_tissue_mapping.get(diagnosis, [])
        
        if not tissues:
            logger.warning(f"No tissue mapping found for diagnosis: {diagnosis}")
        
        logger.info(f"Target tissues for {diagnosis}: {tissues}")
        return tissues

    async def _get_gtex_expression_data_single_queries(self, genes: List[str], gene_to_gencode: Dict[str, str], tissues: List[str]) -> Dict:
        """
        Fetch tissue-specific expression data from GTEx API v2 using single gene-single tissue queries.
        Implements version reduction logic to find working GENCODE IDs.
        """
        gtex_results = {gene: {} for gene in genes}
        
        if not gene_to_gencode or not tissues:
            logger.warning("No valid GENCODE IDs or tissues to query GTEx.")
            return gtex_results

        # Try each gene-tissue combination with version reduction
        for gene, original_gencode_id in gene_to_gencode.items():
            logger.info(f"Querying GTEx for gene {gene} with GENCODE ID {original_gencode_id}")
            
            # Generate version variants
            gencode_variants = self._generate_version_variants(original_gencode_id, max_attempts=10)
            
            # Try each tissue for this gene
            for tissue in tissues:
                logger.info(f"Querying tissue: {tissue}")
                
                # Try each version variant for this gene-tissue combination
                for variant_gencode_id in gencode_variants:
                    logger.debug(f"Trying GENCODE ID variant: {variant_gencode_id} for tissue {tissue}")
                    
                    # Query GTEx with single gene-single tissue format
                    expression_data = await self._fetch_gene_expression_single_gene_tissue(variant_gencode_id, tissue)
                    
                    if expression_data and 'data' in expression_data and expression_data['data']:
                        logger.info(f"Found expression data for {gene} in {tissue} using GENCODE ID: {variant_gencode_id}")
                        
                        # Process the expression data
                        # For single gene-single tissue query, the response structure is different
                        if isinstance(expression_data['data'], list) and len(expression_data['data']) > 0:
                            # Single gene-single tissue returns a list with one record
                            record = expression_data['data'][0]
                            expression_values = record.get('data', [])
                            if expression_values:
                                median_tpm = statistics.median(expression_values)
                                gtex_results[gene][tissue] = {
                                    'median_tpm': median_tpm,
                                    'samples': len(expression_values),
                                    'tissue_found': True,
                                    'tissue_name': tissue,
                                    'gencode_id_used': variant_gencode_id
                                }
                        else:
                            logger.debug(f"Unexpected data structure for {gene} in {tissue}: {expression_data}")
                        
                        # Found data for this tissue, no need to try other variants
                        break
                    else:
                        logger.debug(f"No expression data found for {gene} in {tissue} with GENCODE ID: {variant_gencode_id}")
                else:
                    logger.warning(f"No expression data found for {gene} in {tissue} after trying all version variants")

        logger.info(f"GTEx data retrieval complete for {len(genes)} genes with single gene-tissue queries.")
        return gtex_results

    async def _fetch_gene_expression_single_gene_tissue(self, gencode_id: str, tissue: str) -> Optional[Dict]:
        """
        Fetch gene expression data for a single gene and single tissue from GTEx API v2.
        Uses the specific format: ?gencodeId=X&tissueSiteDetailId=Y
        
        Args:
            gencode_id: GENCODE gene ID (ENSG format)
            tissue: Single tissue name to query
            
        Returns:
            Expression data dictionary or None if failed
        """
        url = f"{self.gtex_api_base}/expression/geneExpression"
        params = {
            'gencodeId': gencode_id,
            'tissueSiteDetailId': tissue,
            'format': 'json'
        }
        
        for attempt in range(self.config['gtex_max_retries']):
            try:
                async with self.session.get(url, params=params) as response: # type: ignore
                    if response.status == 200:
                        content_type = response.headers.get('content-type', '')
                        if 'application/json' in content_type:
                            data = await response.json()
                            logger.debug(f"Successfully fetched expression data for {gencode_id} in {tissue}")
                            return data
                        else:
                            logger.warning(f"GTEx API returned non-JSON response for gene {gencode_id}")
                    elif response.status == 422:
                        error_details = await response.json()
                        logger.debug(f"GTEx API returned 422 for gene {gencode_id}: {error_details}")
                        return None
                    else:
                        logger.warning(f"GTEx API returned status {response.status} for gene {gencode_id}")
                        if attempt < self.config['gtex_max_retries'] - 1:
                            await asyncio.sleep(2 ** attempt)
                        
            except Exception as e:
                logger.warning(f"Error fetching GTEx data for {gencode_id} (attempt {attempt + 1}): {e}")
                if attempt < self.config['gtex_max_retries'] - 1:
                    await asyncio.sleep(2 ** attempt)
        
        logger.error(f"Failed to fetch GTEx data for {gencode_id} after {self.config['gtex_max_retries']} attempts")
        return None

    def _calculate_cross_scale_consistency(self, 
                                         processed_input: Dict,
                                         clinical_data: Dict,
                                         gtex_data: Dict) -> float:
        """
        Calculate coherence across biological scales
        
        Args:
            processed_input: Processed input data
            clinical_data: Clinical data
            gtex_data: GTEx expression data
            
        Returns:
            Cross-scale consistency score (0-1)
        """
        weights = self.config['evidence_weights']
        
        # 1. Variant-level consistency
        variant_score = self._assess_variant_consistency(processed_input)
        
        # 2. Tissue-specific biomarker relevance (GTEx)
        tissue_score = self._assess_tissue_specificity(clinical_data, gtex_data)
        
        # Weighted combination
        consistency_score = (
            weights['variant_consistency'] * variant_score +
            weights['tissue_specificity'] * tissue_score
        )
        
        logger.info(f"Consistency components - Variant: {variant_score:.3f}, Tissue: {tissue_score:.3f}")
        
        return min(1.0, max(0.0, consistency_score))
    
    def _assess_variant_consistency(self, processed_input: Dict) -> float:
        """Assess consistency of variant-level evidence"""
        if 'missense_variants' not in processed_input:
            logger.warning("No missense_variants found in processed input, using default score of 0.5")
            return 0.5
        
        variants = processed_input['missense_variants']
        if not variants:
            logger.warning("Empty missense_variants list, using default score of 0.5")
            return 0.5
        
        pathogenic_scores = []
        clinvar_pathogenic = 0
        alphamissense_pathogenic = 0
        
        for variant in variants:
            # Pathogenicity score
            path_score = variant.get('pathogenicity_score')
            if path_score is not None:
                pathogenic_scores.append(path_score)
            
            # ClinVar annotation
            clinvar_ann = variant.get('clinvar_annotation', '').lower()
            if 'pathogenic' in clinvar_ann:
                clinvar_pathogenic += 1
            
            # AlphaMissense annotation
            alphamissense_ann = variant.get('alphamissense_annotation', '').lower()
            if 'pathogenic' in alphamissense_ann:
                alphamissense_pathogenic += 1
        
        # Calculate consistency
        score = 0.0
        
        # Pathogenicity scores consistency
        if pathogenic_scores:
            avg_pathogenicity = statistics.mean(pathogenic_scores)
            score += 0.4 * avg_pathogenicity
        else:
            logger.warning("No pathogenicity scores found in variants")
        
        # ClinVar consistency
        if variants:
            clinvar_consistency = clinvar_pathogenic / len(variants)
            score += 0.3 * clinvar_consistency
        
        # AlphaMissense consistency
        if variants:
            alphamissense_consistency = alphamissense_pathogenic / len(variants)
            score += 0.3 * alphamissense_consistency
        
        return min(1.0, max(0.0, score))
    
    def _assess_tissue_specificity(self, clinical_data: Dict, gtex_data: Dict) -> float:
        """Assess tissue-specific biomarker relevance using GTEx data"""
        if not gtex_data:
            logger.warning("No GTEx data available, using default tissue specificity score of 0.5")
            return 0.5
        
        diagnosis = clinical_data.get('diagnosis', '').lower()
        biomarkers = clinical_data.get('biomarkers', {})
        
        # Get relevant tissues
        target_tissues = self.disease_tissue_mapping.get(diagnosis, [])
        if not target_tissues:
            logger.warning(f"No tissue mapping found for diagnosis '{diagnosis}', using default score of 0.5")
            return 0.5
        
        tissue_scores = []
        
        for gene, tissue_data in gtex_data.items():
            gene_tissue_scores = []
            
            for tissue in target_tissues:
                if tissue in tissue_data:
                    median_tpm = tissue_data[tissue].get('median_tpm', 0.0)
                    samples = tissue_data[tissue].get('samples', 0)
                    
                    # Score based on expression level and sample size
                    if samples > 0:
                        expression_score = min(1.0, median_tpm / 10.0)  # Normalize by 10 TPM
                        sample_weight = min(1.0, samples / 100.0)  # Weight by sample size
                        tissue_score = expression_score * (0.8 + 0.2 * sample_weight)
                        gene_tissue_scores.append(tissue_score)
            
            if gene_tissue_scores:
                # Use maximum tissue score for each gene
                tissue_scores.append(max(gene_tissue_scores))
        
        # Biomarker-specific adjustments
        if biomarkers and diagnosis == 'breast_cancer':
            # Adjust for breast cancer biomarkers
            her2_status = biomarkers.get('HER2', '').lower()
            er_status = biomarkers.get('ER', '').lower()
            
            if 'positive' in her2_status or 'positive' in er_status:
                # Boost score for positive biomarkers
                tissue_scores = [min(1.0, score * 1.2) for score in tissue_scores]
        
        if not tissue_scores:
            logger.warning("No tissue scores calculated, using default score of 0.5")
            return 0.5
            
        return statistics.mean(tissue_scores)
    
    def _calculate_tissue_specificity_scores(self, clinical_data: Dict, gtex_data: Dict) -> Dict[str, float]:
        """Calculate detailed tissue specificity scores"""
        scores = {}
        
        diagnosis = clinical_data.get('diagnosis', '').lower()
        target_tissues = self.disease_tissue_mapping.get(diagnosis, [])
        
        for gene, tissue_data in gtex_data.items():
            gene_scores = {}
            
            for tissue in target_tissues:
                if tissue in tissue_data:
                    median_tpm = tissue_data[tissue].get('median_tpm', 0.0)
                    samples = tissue_data[tissue].get('samples', 0)
                    
                    # Calculate tissue-specific score
                    if samples > 0:
                        score = min(1.0, median_tpm / 10.0)
                        gene_scores[tissue] = score
            
            scores[gene] = gene_scores
        
        return scores
    
    def export_results(self, result: CoherenceResult, output_file: str) -> None:
        """Export coherence results to JSON file"""
        output_data = {
            'cross_scale_consistency': result.cross_scale_consistency,
            'tissue_specificity_scores': result.tissue_specificity_scores,
            'evidence_weights': result.evidence_weights,
            'gtex_validation_results': result.gtex_validation_results,
            'processing_timestamp': result.processing_timestamp
        }
        
        try:
            with open(output_file, 'w') as f:
                json.dump(output_data, f, indent=2)
            logger.info(f"Coherence results exported to: {output_file}")
        except Exception as e:
            logger.error(f"Failed to export results to {output_file}: {e}")
            raise



async def main():
    """Main function for command-line usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Coherence Analyzer")
    parser.add_argument("--processed-input", required=True,
                       help="Processed input file (required)")
    parser.add_argument("--clinical-data", required=True,
                       help="Clinical data file (required)")
    parser.add_argument("-o", "--output", default="coherence_results.json",
                       help="Output file (default: coherence_results.json)")
    parser.add_argument("--simple", action="store_true",
                       help="Output only cross_scale_consistency score")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Warn about default output file
    if args.output == "coherence_results.json":
        logger.warning("Using default output file: coherence_results.json")
    
    # Initialize analyzer
    analyzer = CoherenceAnalyzer()
    
    # Run analysis with provided files
    result = await analyzer.analyze_coherence(args.processed_input, args.clinical_data)
    
    if args.simple:
        # Simple output format as requested
        simple_output = {"cross_scale_consistency": result.cross_scale_consistency}
        logger.info(json.dumps(simple_output, indent=2))
    else:
        # Full output
        analyzer.export_results(result, args.output)
        logger.info(f"✓ Coherence analysis complete")
        logger.info(f"✓ Results saved to: {args.output}")
        logger.info(f"✓ Cross-scale consistency: {result.cross_scale_consistency:.3f}")


if __name__ == "__main__":
    asyncio.run(main()) 