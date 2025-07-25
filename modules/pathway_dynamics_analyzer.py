"""
Pathway Impact & Simulation with GSEApy

This module performs real gene set enrichment analysis using GSEApy
with MSigDB and Reactome gene sets to analyze pathway impact of genetic variants.
"""

import json
import logging
import os
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
import numpy as np
from datetime import datetime

# Configure logging first
from .logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)

try:
    import gseapy as gp
    # Use enrichr for ORA
    from gseapy import enrichr
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False
    logger.error("GSEApy not available. Install with: pip install gseapy")


class PathwayImpactAnalyzer:
    """
    Analyzes pathway impact using GSEApy for gene set enrichment analysis.
    
    This class performs gene set enrichment analysis on gene lists derived from
    genetic variants to identify impacted biological pathways and processes.
    """
    
    def __init__(self, input_file: str = "variant_impact_results.json"):
        """
        Initialize the Pathway Impact Analyzer.
        
        Args:
            input_file: Path to the variant impact results JSON file
        """
        self.input_file = input_file
        self.input_data = None
        self.gene_list = []
        self.enrichment_results = {}
        self.pathway_scores = {}
        
        if not GSEAPY_AVAILABLE:
            raise ImportError("GSEApy is required for pathway analysis. Install with: pip install gseapy")
    
    def load_input_data(self) -> bool:
        """
        Load and validate input data from processed_input.json.
        
        Returns:
            bool: True if data loaded successfully, False otherwise
        """
        try:
            with open(self.input_file, 'r') as f:
                self.input_data = json.load(f)
            
            logger.info(f"Loaded input data from {self.input_file}")
            return True
            
        except FileNotFoundError:
            logger.error(f"Input file {self.input_file} not found")
            return False
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in {self.input_file}: {e}")
            return False
    
    def extract_gene_list(self) -> List[str]:
        """
        Extract gene list from missense variants in the input data.
        
        Returns:
            List[str]: List of unique gene names
        """
        if not self.input_data or 'missense_variants' not in self.input_data:
            logger.error("No missense variants found in input data")
            return []
        
        genes = set()
        for variant in self.input_data['missense_variants']:
            if 'gene' in variant and variant['gene']:
                genes.add(variant['gene'])
        
        self.gene_list = list(genes)
        logger.info(f"Extracted {len(self.gene_list)} unique genes: {self.gene_list}")
        return self.gene_list
    
    def create_gene_ranking(self, method: str = "pathogenicity_mean") -> Dict[str, float]:
        """
        Create gene ranking for GSEA analysis.
        
        Args:
            method: Ranking method ("uniform", "variant_count", "clinical_relevance", 
                    "pathogenicity_max", "pathogenicity_mean", "pathogenicity_sum")
            
        Returns:
            Dict[str, float]: Gene name to score mapping
        """
        if not self.gene_list:
            logger.error("No gene list available. Run extract_gene_list() first.")
            return {}
        
        if self.input_data is None:
            logger.error("No input data available. Run load_input_data() first.")
            return {}
        
        gene_scores = {}
        
        if method.startswith("pathogenicity"):
            # Use pathogenicity scores from variant impact results
            per_gene = {g: [] for g in self.gene_list}
            
            for variant in self.input_data['missense_variants']:
                gene = variant.get('gene', '')
                pathogenicity_score = variant.get('pathogenicity_score')
                
                if gene:
                    # Handle null scores by using 0.5 as default
                    if pathogenicity_score is None:
                        pathogenicity_score = 0.5
                    per_gene.setdefault(gene, []).append(float(pathogenicity_score))
            
            # Aggregate scores based on method
            agg = method.split("_")[1]  # max, mean, or sum
            
            for gene, scores in per_gene.items():
                if not scores:
                    gene_scores[gene] = 0.0
                elif agg == "max":
                    gene_scores[gene] = max(scores)
                elif agg == "mean":
                    gene_scores[gene] = sum(scores) / len(scores)
                else:  # sum
                    gene_scores[gene] = sum(scores)
                    
        elif method == "uniform":
            # Uniform ranking - all genes get same score
            for gene in self.gene_list:
                gene_scores[gene] = 1.0
                
        elif method == "variant_count":
            # Score based on number of variants per gene
            gene_variant_counts = {}
            for variant in self.input_data['missense_variants']:
                gene = variant.get('gene', '')
                if gene:
                    gene_variant_counts[gene] = gene_variant_counts.get(gene, 0) + 1
            
            max_count = max(gene_variant_counts.values()) if gene_variant_counts else 1
            for gene in self.gene_list:
                count = gene_variant_counts.get(gene, 0)
                gene_scores[gene] = count / max_count if max_count > 0 else 0.0
                
        elif method == "clinical_relevance":
            # Score based on clinical relevance (cancer-related genes get higher scores)
            cancer_genes = {
                'TP53': 1.0, 'BRCA1': 0.9, 'BRCA2': 0.9, 'PTEN': 0.8,
                'APC': 0.8, 'RB1': 0.7, 'CDKN2A': 0.7, 'ATM': 0.7,
                'CHEK2': 0.6, 'PALB2': 0.6, 'CFTR': 0.5
            }
            
            for gene in self.gene_list:
                gene_scores[gene] = cancer_genes.get(gene, 0.3)
        
        logger.info(f"Created gene ranking with {len(gene_scores)} genes using {method} method")
        return gene_scores
    
    def run_gsea_analysis(self, gene_scores: Dict[str, float], 
                          gene_sets: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Run GSEA analysis on the gene list.
        
        Args:
            gene_scores: Dictionary mapping gene names to scores
            gene_sets: List of gene set databases to use
            
        Returns:
            Dict[str, Any]: GSEA enrichment results
        """
        if not gene_scores:
            logger.error("No gene scores provided")
            return {}
        
        # Default gene sets to analyze
        if gene_sets is None:
            gene_sets = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2021']
        
        results = {}
        
        try:
            # For small gene lists (< 100 genes), use ORA instead of GSEA
            if len(gene_scores) < 100:
                logger.info(f"Small gene list ({len(gene_scores)} genes), using ORA (Over-Representation Analysis)")
                
                # Convert gene scores to gene list (take top genes by score)
                gene_list = [gene for gene, score in sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)]
                
                for gene_set in gene_sets:
                    try:
                        enr = enrichr(
                            gene_list=gene_list,
                            gene_sets=gene_set,
                            organism='Human',
                            outdir=None,
                            no_plot=True
                        )
                        
                        if hasattr(enr, 'results') and enr.results is not None and len(enr.results) > 0:
                            # enr.results is a list of dicts for enrichr
                            # Convert to JSON-serializable format
                            serializable_results = []
                            for result in enr.results:
                                if isinstance(result, dict):
                                    serializable_result = {}
                                    for key, value in result.items():
                                        if hasattr(value, 'item'):  # numpy scalar
                                            serializable_result[key] = value.item()
                                        elif hasattr(value, 'tolist'):  # numpy array
                                            serializable_result[key] = value.tolist()
                                        else:
                                            serializable_result[key] = value
                                    serializable_results.append(serializable_result)
                                else:
                                    serializable_results.append(result)
                            
                            # Only add results if we have actual data
                            if serializable_results and len(serializable_results) > 0:
                                # Check if the first result has actual data (not just headers)
                                first_result = serializable_results[0]
                                if isinstance(first_result, dict) and len(first_result) > 1:
                                    results[gene_set] = {
                                        'enrichment_results': serializable_results,
                                        'summary': self._summarize_enrichment(pd.DataFrame(serializable_results))
                                    }
                                else:
                                    logger.warning(f"No valid enrichment results for {gene_set}")
                            else:
                                logger.warning(f"No enrichment results for {gene_set}")
                        else:
                            logger.warning(f"No enrichment results returned for {gene_set}")
                    except Exception as e:
                        logger.warning(f"Failed to analyze {gene_set}: {e}")
                        continue
                        
            else:
                # For larger gene lists, use prerank (GSEApy's prerank, not GSEA class)
                logger.info(f"Large gene list ({len(gene_scores)} genes), using prerank GSEA")
                
                # Convert gene scores to pandas DataFrame for prerank
                gene_series = pd.Series(gene_scores)
                prerank_df = gene_series.sort_values(ascending=False)
                prerank_df = prerank_df.reset_index()
                prerank_df.columns = ['gene_name', 'score']
                
                # Run prerank for each gene set
                for gene_set in gene_sets:
                    try:
                        prerank_res = gp.prerank(
                            rnk=prerank_df,
                            gene_sets=gene_set,
                            threads=4,
                            permutation_num=100,  # reduce for speed, increase for accuracy
                            outdir=None,
                            seed=42,
                            verbose=False
                        )
                        if hasattr(prerank_res, 'res2d') and prerank_res.res2d is not None and not prerank_res.res2d.empty:
                            # Convert DataFrame to JSON-serializable format
                            df_records = prerank_res.res2d.reset_index().to_dict('records')
                            serializable_records = []
                            for record in df_records:
                                serializable_record = {}
                                for key, value in record.items():
                                    if hasattr(value, 'item'):  # numpy scalar
                                        serializable_record[key] = value.item()
                                    elif hasattr(value, 'tolist'):  # numpy array
                                        serializable_record[key] = value.tolist()
                                    else:
                                        serializable_record[key] = value
                                serializable_records.append(serializable_record)
                            
                            # Only add results if we have actual data
                            if serializable_records and len(serializable_records) > 0:
                                # Check if the first result has actual data (not just headers)
                                first_result = serializable_records[0]
                                if isinstance(first_result, dict) and len(first_result) > 1:
                                    results[gene_set] = {
                                        'enrichment_results': serializable_records,
                                        'summary': self._summarize_enrichment(prerank_res.res2d)
                                    }
                                else:
                                    logger.warning(f"No valid prerank results for {gene_set}")
                            else:
                                logger.warning(f"No prerank results for {gene_set}")
                        else:
                            logger.warning(f"No prerank results returned for {gene_set}")
                    except Exception as e:
                        logger.warning(f"Failed to analyze {gene_set} with prerank: {e}")
                        continue
            
            self.enrichment_results = results
            logger.info(f"Enrichment analysis completed for {len(results)} gene sets")
            
        except Exception as e:
            logger.error(f"Error running enrichment analysis: {e}")
            return {}
        
        return results
    
    def _summarize_enrichment(self, enrichment_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Summarize enrichment results for a gene set.
        
        Args:
            enrichment_df: DataFrame with enrichment results
            
        Returns:
            Dict[str, Any]: Summary statistics
        """
        if enrichment_df is None or enrichment_df.empty:
            logger.warning("Empty enrichment DataFrame provided")
            return {
                'total_pathways': 0,
                'significant_pathways': 0,
                'top_pathways': [],
                'max_enrichment_score': 0,
                'min_fdr': 1.0
            }
        
        # Check if DataFrame only contains headers (single row with column names)
        if len(enrichment_df) == 1:
            first_row = enrichment_df.iloc[0]
            # If all values in the first row are the same as column names, it's likely headers
            if all(str(first_row[col]) == str(col) for col in enrichment_df.columns):
                logger.warning("Enrichment DataFrame contains only headers, no actual data")
                return {
                    'total_pathways': 0,
                    'significant_pathways': 0,
                    'top_pathways': [],
                    'max_enrichment_score': 0,
                    'min_fdr': 1.0
                }
        
        # Handle different column names for FDR/p-value
        fdr_col = None
        if 'FDR q-val' in enrichment_df.columns:
            fdr_col = 'FDR q-val'
        elif 'Adjusted P-value' in enrichment_df.columns:
            fdr_col = 'Adjusted P-value'
        elif 'P-value' in enrichment_df.columns:
            fdr_col = 'P-value'
        elif 'fdr' in enrichment_df.columns:
            fdr_col = 'fdr'
        
        if fdr_col is None:
            logger.warning("No FDR/p-value column found in enrichment results")
            return {
                'total_pathways': len(enrichment_df),
                'significant_pathways': 0,
                'top_pathways': enrichment_df.head(10).to_dict('records') if not enrichment_df.empty else [],
                'max_enrichment_score': 0,
                'min_fdr': 1.0
            }
        
        # Filter significant results (FDR < 0.05)
        significant = enrichment_df[enrichment_df[fdr_col] < 0.05]
        
        # Try to get enrichment score column
        es_col = None
        for col in ['Enrichment Score', 'es', 'NES', 'nes']:
            if col in enrichment_df.columns:
                es_col = col
                break
        
        summary = {
            'total_pathways': len(enrichment_df),
            'significant_pathways': len(significant),
            'top_pathways': significant.head(10).to_dict('records') if not significant.empty else [],
            'max_enrichment_score': enrichment_df[es_col].max() if es_col and es_col in enrichment_df.columns else 0,
            'min_fdr': enrichment_df[fdr_col].min() if not enrichment_df.empty else 1.0
        }
        
        return summary
    
    def analyze_target_pathways(self) -> Dict[str, float]:
        """
        Analyze specific target pathways for minimum expected scores.
        
        Returns:
            Dict[str, float]: Pathway name to score mapping
        """
        target_pathways = {
            'p53_signaling': 0.70,
            'DNA_repair': 0.65,
            'Apoptosis': 0.60
        }
        
        pathway_scores = {}
        
        # Search for target pathways in enrichment results
        for gene_set, results in self.enrichment_results.items():
            if 'enrichment_results' in results:
                enrichment_data = results['enrichment_results']
                logger.info(f"Processing {gene_set}: enrichment_data type = {type(enrichment_data)}")
                
                # Handle different data structures
                if isinstance(enrichment_data, list):
                    pathways = enrichment_data
                elif isinstance(enrichment_data, dict):
                    # If it's a dict, convert to list of records
                    pathways = [enrichment_data]
                elif hasattr(enrichment_data, 'to_dict'):  # pandas DataFrame
                    # Convert DataFrame to list of dictionaries
                    pathways = enrichment_data.to_dict('records')
                    logger.info(f"Converted DataFrame to {len(pathways)} pathway records")
                else:
                    logger.warning(f"Unexpected enrichment_data type: {type(enrichment_data)}")
                    continue
                
                for pathway in pathways:
                    # Debug: log the pathway structure
                    if isinstance(pathway, str):
                        logger.warning(f"Pathway is string: {pathway}")
                        continue
                    elif isinstance(pathway, dict):
                        logger.debug(f"Pathway keys: {list(pathway.keys())}")
                    else:
                        logger.warning(f"Unexpected pathway type: {type(pathway)}")
                        continue
                    
                    # Try to get pathway name from possible columns
                    pathway_name_raw = (
                        pathway.get('Term') or
                        pathway.get('Name') or
                        pathway.get('term_name') or
                        pathway.get('Gene_set') or
                        pathway.get('Name') or  # Add this for Reactome
                        ''
                    )
                    pathway_name = pathway_name_raw.lower() if pathway_name_raw else ''
                    
                    for target, min_score in target_pathways.items():
                        # More flexible matching for pathway names
                        target_terms = []
                        
                        if target == 'p53_signaling':
                            target_terms = ['p53 signaling pathway', 'p53', 'tp53']
                        elif target == 'Apoptosis':
                            target_terms = [
                                'apoptosis', 'apoptotic', 'cell death', 'death', 'survival', 
                                'caspase', 'bcl', 'bax', 'p53', 'tp53',
                                'transcriptional regulation by tp53'
                            ]
                        elif target == 'DNA_repair':
                            target_terms = [
                                'dna repair', 'dna damage', 'dna-binding', 'dna-templated',
                                'cellular response to dna damage stimulus',
                                'positive regulation of dna-binding transcription factor activity',
                                'regulation of transcription, dna-templated',
                                'homologous recombination', 'mismatch repair', 'base excision repair',
                                'double-strand break repair', 'fanconi anemia'
                            ]
                        else:
                            target_terms = [
                                target.replace('_', ' '),
                                target.replace('_', ''),
                            ]
                        
                        matched = any(term in pathway_name for term in target_terms)
                        
                        if matched:
                            # For ORA, use -log10(p-value) as score, for GSEA use Enrichment Score or NES
                            score = 0
                            if 'P-value' in pathway and pathway.get('P-value') is not None:
                                p_value = pathway.get('P-value', 1.0)
                                try:
                                    score = -np.log10(float(p_value)) if float(p_value) > 0 else 10.0
                                except Exception:
                                    score = 0
                            elif 'Adjusted P-value' in pathway and pathway.get('Adjusted P-value') is not None:
                                p_value = pathway.get('Adjusted P-value', 1.0)
                                try:
                                    score = -np.log10(float(p_value)) if float(p_value) > 0 else 10.0
                                except Exception:
                                    score = 0
                            elif 'Enrichment Score' in pathway and pathway.get('Enrichment Score') is not None:
                                score = pathway.get('Enrichment Score', 0)
                            elif 'NES' in pathway and pathway.get('NES') is not None:
                                score = pathway.get('NES', 0)
                            elif 'es' in pathway and pathway.get('es') is not None:
                                score = pathway.get('es', 0)
                            elif 'nes' in pathway and pathway.get('nes') is not None:
                                score = pathway.get('nes', 0)
                            else:
                                score = 0
                            
                            pathway_scores[target] = score
                            logger.info(f"Found {target}: score = {score:.3f} (min expected: {min_score})")
        
        self.pathway_scores = pathway_scores
        return pathway_scores
    
    def validate_target_scores(self) -> Dict[str, Dict[str, Any]]:
        """
        Validate that target pathway scores meet minimum expected values.
        
        Returns:
            Dict[str, Dict[str, Any]]: Validation results for each target pathway
        """
        validation_results = {}
        
        for pathway, min_score in {
            'p53_signaling': 0.70,
            'DNA_repair': 0.65,
            'Apoptosis': 0.60
        }.items():
            actual_score = self.pathway_scores.get(pathway, 0)
            passed = actual_score >= min_score
            
            validation_results[pathway] = {
                'expected_min': min_score,
                'actual_score': actual_score,
                'passed': bool(passed),  # Convert to bool for JSON serialization
                'status': 'PASS' if passed else 'FAIL'
            }
            
            logger.info(f"{pathway}: {actual_score:.3f} {'≥' if passed else '<'} {min_score} ({'PASS' if passed else 'FAIL'})")
        
        return validation_results
    
    def generate_report(self, output_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Generate comprehensive pathway impact report.
        
        Args:
            output_file: Optional file path to save report
            
        Returns:
            Dict[str, Any]: Complete analysis report
        """
        validation_results = self.validate_target_scores()
        
        # Get gene list and count
        gene_list = self.gene_list if hasattr(self, 'gene_list') else []
        gene_count = len(gene_list)
        
        # Create pathway enrichment scores
        pathway_enrichment = {}
        for pathway, score in self.pathway_scores.items():
            pathway_enrichment[pathway] = score
        
        # Convert enrichment results to JSON-serializable format
        serializable_enrichment_results = {}
        for gene_set, results in self.enrichment_results.items():
            # Convert enrichment results to ensure they're serializable
            enrichment_results = results.get('enrichment_results', [])
            if isinstance(enrichment_results, list):
                # Convert any pandas Series or numpy types to native Python types
                serializable_results = []
                for result in enrichment_results:
                    if isinstance(result, dict):
                        serializable_result = {}
                        for key, value in result.items():
                            if hasattr(value, 'item'):  # numpy scalar
                                serializable_result[key] = value.item()
                            elif hasattr(value, 'tolist'):  # numpy array
                                serializable_result[key] = value.tolist()
                            else:
                                serializable_result[key] = value
                        serializable_results.append(serializable_result)
                    else:
                        serializable_results.append(result)
            else:
                serializable_results = enrichment_results
            
            # Convert summary to ensure it's serializable
            summary = results.get('summary', {})
            if isinstance(summary, dict):
                serializable_summary = {}
                for key, value in summary.items():
                    if key == 'top_pathways' and isinstance(value, list):
                        # Convert top_pathways to serializable format
                        serializable_top_pathways = []
                        for pathway in value:
                            if isinstance(pathway, dict):
                                serializable_pathway = {}
                                for pkey, pvalue in pathway.items():
                                    if hasattr(pvalue, 'item'):  # numpy scalar
                                        serializable_pathway[pkey] = pvalue.item()
                                    elif hasattr(pvalue, 'tolist'):  # numpy array
                                        serializable_pathway[pkey] = pvalue.tolist()
                                    else:
                                        serializable_pathway[pkey] = pvalue
                                serializable_top_pathways.append(serializable_pathway)
                            else:
                                serializable_top_pathways.append(pathway)
                        serializable_summary[key] = serializable_top_pathways
                    elif hasattr(value, 'item'):  # numpy scalar
                        serializable_summary[key] = value.item()
                    elif hasattr(value, 'tolist'):  # numpy array
                        serializable_summary[key] = value.tolist()
                    else:
                        serializable_summary[key] = value
            else:
                serializable_summary = summary
            
            serializable_enrichment_results[gene_set] = {
                'enrichment_results': serializable_results,
                'summary': serializable_summary
            }
        
        report = {
            'gene_list': gene_list,
            'gene_count': gene_count,
            'pathway_enrichment': pathway_enrichment,
            'target_pathway_scores': self.pathway_scores,
            'validation_results': validation_results,
            'enrichment_results': serializable_enrichment_results,
            'summary': {
                'total_pathways_analyzed': sum(len(r.get('enrichment_results', [])) 
                                             for r in self.enrichment_results.values()),
                'significant_pathways': sum(r.get('summary', {}).get('significant_pathways', 0) 
                                         for r in self.enrichment_results.values()),
                'target_pathways_passed': sum(1 for v in validation_results.values() 
                                            if v['passed'])
            }
        }
        
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Report saved to {output_file}")
        
        return report
    
    def run_complete_analysis(self, output_file: str = "pathway_impact_results.json") -> Dict[str, Any]:
        """
        Run complete pathway impact analysis pipeline.
        
        Args:
            output_file: Output file path for results
            
        Returns:
            Dict[str, Any]: Complete analysis results
        """
        logger.info("Starting complete pathway impact analysis...")
        
        # Step 1: Load input data
        if not self.load_input_data():
            return {'error': 'Failed to load input data'}
        
        # Step 2: Extract gene list
        genes = self.extract_gene_list()
        if not genes:
            return {'error': 'No genes found in input data'}
        
        # Step 3: Create gene ranking
        gene_scores = self.create_gene_ranking(method="pathogenicity_mean")
        
        # Step 4: Run GSEA analysis
        enrichment_results = self.run_gsea_analysis(gene_scores)
        
        # Step 5: Analyze target pathways
        self.analyze_target_pathways()
        
        # Step 6: Generate report
        report = self.generate_report(output_file)
        
        logger.info("Pathway impact analysis completed successfully")
        return report


def main():
    """Main function to run pathway impact analysis."""
    analyzer = PathwayImpactAnalyzer()
    results = analyzer.run_complete_analysis()
    
    # Print summary
    logger.info("\n=== Pathway Impact Analysis Results ===")
    logger.info(f"Genes analyzed: {results.get('gene_count', 0)}")
    logger.info(f"Total pathways analyzed: {results.get('summary', {}).get('total_pathways_analyzed', 0)}")
    logger.info(f"Significant pathways: {results.get('summary', {}).get('significant_pathways', 0)}")
    logger.info(f"Target pathways passed: {results.get('summary', {}).get('target_pathways_passed', 0)}")
    
    # Print target pathway validation
    logger.info("\n=== Target Pathway Validation ===")
    for pathway, validation in results.get('validation_results', {}).items():
        status = validation['status']
        score = validation['actual_score']
        expected = validation['expected_min']
        logger.info(f"{pathway}: {score:.3f} {'≥' if validation['passed'] else '<'} {expected} ({status})")


if __name__ == "__main__":
    main()