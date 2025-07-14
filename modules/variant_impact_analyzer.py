#!/usr/bin/env python3
"""
Variant Impact Analyzer

A module for analyzing variant impact using AlphaMissense and ClinVar databases.
Performs parallel processing of variants and their annotations.

Example usage:
    analyzer = VariantImpactAnalyzer()
    results = analyzer.analyze_variants("processed_input.json")
"""

import os
import json
import sqlite3
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import statistics
from collections import Counter

# Import the ClinVar annotator
try:
    from .clinvar_annotator import ClinVarAnnotator, ClinVarAnnotation
except ImportError:
    from clinvar_annotator import ClinVarAnnotator, ClinVarAnnotation

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class AlphaMissenseResult:
    """AlphaMissense annotation result"""
    variant_id: str
    gene: str
    protein_change: str
    average_pathogenicity_score: Optional[float] = None
    max_occurring_class: Optional[str] = None
    matching_transcripts: List[str] = field(default_factory=list)
    total_matches: int = 0
    warnings: List[str] = field(default_factory=list)


@dataclass
class VariantImpactResult:
    """Complete variant impact analysis result"""
    variant_id: str
    gene: str
    protein_change: str
    alphamissense: AlphaMissenseResult
    clinvar: ClinVarAnnotation
    processing_timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class VariantImpactAnalyzer:
    """Analyzer for variant impact using AlphaMissense and ClinVar"""
    
    def __init__(self, alphamissense_db_path: str = "cache/alphamissense/alphamissense_hg38.db"):
        """
        Initialize the variant impact analyzer
        
        Args:
            alphamissense_db_path: Path to AlphaMissense database
        """
        self.alphamissense_db_path = Path(alphamissense_db_path)
        self.clinvar_annotator = ClinVarAnnotator()
        
        # Validate database exists
        if not self.alphamissense_db_path.exists():
            raise FileNotFoundError(f"AlphaMissense database not found: {self.alphamissense_db_path}")
        
        logger.info("Variant Impact Analyzer initialized")
    
    def analyze_variants(self, input_file: str) -> Dict[str, Any]:
        """
        Analyze all variants in the processed input file
        
        Args:
            input_file: Path to processed_input.json file
            
        Returns:
            Dictionary containing analysis results and original data
        """
        # Load input data
        with open(input_file, 'r') as f:
            input_data = json.load(f)
        
        variants = input_data.get('missense_variants', [])
        logger.info(f"Analyzing {len(variants)} variants")
        
        # Process each variant and add results directly to the variant
        for variant in variants:
            variant_id = variant['id']
            logger.info(f"Processing variant: {variant_id}")
            
            # Analyze variant
            result = self._analyze_single_variant(variant)
            
            # Add AlphaMissense results to variant
            variant['pathogenicity_score'] = result.alphamissense.average_pathogenicity_score
            variant['alphamissense_annotation'] = result.alphamissense.max_occurring_class.title() if result.alphamissense.max_occurring_class else None
            
            # Add ClinVar results to variant
            variant['clinvar_variation_id'] = result.clinvar.variation_id
            variant['clinvar_annotation'] = result.clinvar.clinical_significance.title() if result.clinvar.clinical_significance else None
        
        # Add analysis metadata
        input_data['variant_impact_analysis'] = {
            'analysis_timestamp': datetime.now().isoformat(),
            'total_variants': len(variants)
        }
        
        return input_data
    
    def _analyze_single_variant(self, variant: Dict) -> VariantImpactResult:
        """
        Analyze a single variant using AlphaMissense and ClinVar
        
        Args:
            variant: Variant dictionary from input
            
        Returns:
            VariantImpactResult with analysis results
        """
        variant_id = variant['id']
        gene = variant['gene']
        protein_change = variant['protein_change']
        vep_transcript_ids = variant.get('vep_transcript_ids', [])
        
        # Task 1: AlphaMissense analysis
        alphamissense_result = self._analyze_alphamissense(
            variant_id, gene, protein_change, vep_transcript_ids
        )
        
        # Task 2: ClinVar analysis
        clinvar_result = self.clinvar_annotator.annotate_variant(variant)
        
        return VariantImpactResult(
            variant_id=variant_id,
            gene=gene,
            protein_change=protein_change,
            alphamissense=alphamissense_result,
            clinvar=clinvar_result
        )
    
    def _analyze_alphamissense(self, variant_id: str, gene: str, 
                              protein_change: str, vep_transcript_ids: List[str]) -> AlphaMissenseResult:
        """
        Analyze variant using AlphaMissense database
        
        Args:
            variant_id: Variant identifier
            gene: Gene name
            protein_change: Protein change (e.g., "p.Arg175His")
            vep_transcript_ids: List of transcript IDs to search
            
        Returns:
            AlphaMissenseResult with analysis results
        """
        result = AlphaMissenseResult(
            variant_id=variant_id,
            gene=gene,
            protein_change=protein_change
        )
        
        try:
            # Extract position number from protein change
            position_match = re.search(r'(\d+)', protein_change)
            if not position_match:
                result.warnings.append(f"Could not extract position from protein change: {protein_change}")
                return result
            
            position_number = position_match.group(1)
            
            # Connect to database
            conn = sqlite3.connect(self.alphamissense_db_path)
            cursor = conn.cursor()
            
            # Search for variants in the specified transcripts with matching position
            query = """
                SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant
                FROM alphamissense 
                WHERE transcript_id IN ({})
                AND protein_variant LIKE ?
                ORDER BY am_pathogenicity DESC
            """.format(','.join(['?' for _ in vep_transcript_ids]))
            
            # Create pattern for position matching
            position_pattern = f"%{position_number}%"
            
            # Execute query
            params = vep_transcript_ids + [position_pattern]
            cursor.execute(query, params)
            matches = cursor.fetchall()
            
            # If too many matches, try to filter by exact protein_variant
            if matches and len(matches) > 10:
                # Try to match the full protein_change string (e.g., 'p.Glu23fs')
                strict_query = """
                    SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant
                    FROM alphamissense 
                    WHERE transcript_id IN ({})
                    AND protein_variant = ?
                    ORDER BY am_pathogenicity DESC
                """.format(','.join(['?' for _ in vep_transcript_ids]))
                strict_params = vep_transcript_ids + [protein_change]
                cursor.execute(strict_query, strict_params)
                strict_matches = cursor.fetchall()
                if strict_matches:
                    matches = strict_matches
                else:
                    # If strict matching fails, try filtering by pathogenicity score
                    pathogenicity_query = """
                        SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant
                        FROM alphamissense 
                        WHERE transcript_id IN ({})
                        AND protein_variant LIKE ?
                        AND am_pathogenicity > 0.8
                        ORDER BY am_pathogenicity DESC
                        LIMIT 10
                    """.format(','.join(['?' for _ in vep_transcript_ids]))
                    pathogenicity_params = vep_transcript_ids + [position_pattern]
                    print(pathogenicity_query, pathogenicity_params)
                    cursor.execute(pathogenicity_query, pathogenicity_params)
                    pathogenicity_matches = cursor.fetchall()
                    if pathogenicity_matches:
                        matches = pathogenicity_matches
                        result.warnings.append(f"Filtered to {len(matches)} high-pathogenicity matches (score > 0.8)")
                    else:
                        result.warnings.append(f"Strict filtering by protein_variant and pathogenicity gave no results; using broader position-based matches.")
            
            if matches:
                # Extract pathogenicity scores and classes
                pathogenicity_scores = [match[0] for match in matches if match[0] is not None]
                classes = [match[1] for match in matches if match[1] is not None]
                matching_transcripts = list(set([match[3] for match in matches]))
                
                if pathogenicity_scores:
                    result.average_pathogenicity_score = statistics.mean(pathogenicity_scores)
                
                if classes:
                    # Find most common class
                    class_counter = Counter(classes)
                    result.max_occurring_class = class_counter.most_common(1)[0][0]
                
                result.matching_transcripts = matching_transcripts
                result.total_matches = len(matches)
                
                logger.info(f"Found {len(matches)} AlphaMissense matches for {variant_id}")
            else:
                result.warnings.append(f"No AlphaMissense matches found for {variant_id}")
                logger.warning(f"No AlphaMissense matches found for {variant_id}")
            
            conn.close()
        
        except Exception as e:
            error_msg = f"Error analyzing AlphaMissense for {variant_id}: {str(e)}"
            result.warnings.append(error_msg)
            logger.error(error_msg)
        
        return result
    
    def save_results(self, results: Dict[str, Any], output_file: str):
        """
        Save analysis results to file
        
        Args:
            results: Analysis results dictionary
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        logger.info(f"Results saved to: {output_file}")


def main():
    """Main function for testing"""
    analyzer = VariantImpactAnalyzer()
    
    # Analyze variants
    results = analyzer.analyze_variants("processed_input.json")
    
    # Save results
    analyzer.save_results(results, "variant_impact_results.json")
    
    # Print summary
    print(f"Analysis complete.")
    
    # Print sample results
    for variant_id, result in results.get('variant_impact_analysis', {}).get('results', {}).items():
        print(f"\nVariant: {variant_id}")
        print(f"  Gene: {result.gene}")
        print(f"  Protein Change: {result.protein_change}")
        print(f"  AlphaMissense Score: {result.alphamissense.average_pathogenicity_score}")
        print(f"  AlphaMissense Class: {result.alphamissense.max_occurring_class}")
        print(f"  ClinVar Significance: {result.clinvar.clinical_significance}")
        if result.alphamissense.warnings:
            print(f"  Warnings: {result.alphamissense.warnings}")


if __name__ == "__main__":
    main() 