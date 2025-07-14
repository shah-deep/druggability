#!/usr/bin/env python3
"""
Test script for Pathway Impact Analysis using GSEApy

This script tests the pathway impact analyzer with the large variant file
and validates the minimum expected scores for target pathways.
"""

import sys
import os
import json
import logging
import numpy as np
from pathlib import Path

# Add modules to path
sys.path.append(str(Path(__file__).parent / "modules"))

from pathway_impact_analyzer import PathwayImpactAnalyzer

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def test_pathway_analysis():
    """Test pathway impact analysis with large variant file"""
    
    # Input file
    input_file = "variant_impact_results_large.json"
    
    if not os.path.exists(input_file):
        logger.error(f"Input file {input_file} not found")
        return False
    
    logger.info(f"Starting pathway impact analysis with {input_file}")
    
    try:
        # Initialize analyzer
        analyzer = PathwayImpactAnalyzer(input_file=input_file)
        
        # Run complete analysis
        results = analyzer.run_complete_analysis(output_file="pathway_impact_results.json")
        
        if 'error' in results:
            logger.error(f"Analysis failed: {results['error']}")
            return False
        
        # Print detailed results
        print("\n" + "="*60)
        print("PATHWAY IMPACT ANALYSIS RESULTS")
        print("="*60)
        
        print(f"\nInput file: {results.get('input_file', 'N/A')}")
        print(f"Analysis timestamp: {results.get('analysis_timestamp', 'N/A')}")
        print(f"Genes analyzed: {results.get('gene_count', 0)}")
        print(f"Gene list: {', '.join(results.get('gene_list', []))}")
        
        # Print enrichment summary
        summary = results.get('summary', {})
        print(f"\nEnrichment Summary:")
        print(f"  Total pathways analyzed: {summary.get('total_pathways_analyzed', 0)}")
        print(f"  Significant pathways: {summary.get('significant_pathways', 0)}")
        print(f"  Target pathways passed: {summary.get('target_pathways_passed', 0)}")
        
        # Print target pathway validation
        print(f"\nTarget Pathway Validation:")
        print("-" * 40)
        validation_results = results.get('validation_results', {})
        
        for pathway, validation in validation_results.items():
            status = validation['status']
            score = validation['actual_score']
            expected = validation['expected_min']
            passed = validation['passed']
            
            print(f"{pathway:15} | {score:8.3f} | {'≥' if passed else '<':2} {expected:6.2f} | {status}")
        
        # Print detailed enrichment results
        print(f"\nDetailed Enrichment Results:")
        print("-" * 40)
        
        enrichment_results = results.get('enrichment_results', {})
        for gene_set, results_data in enrichment_results.items():
            print(f"\n{gene_set}:")
            summary_data = results_data.get('summary', {})
            print(f"  Total pathways: {summary_data.get('total_pathways', 0)}")
            print(f"  Significant pathways: {summary_data.get('significant_pathways', 0)}")
            print(f"  Max enrichment score: {summary_data.get('max_enrichment_score', 0):.3f}")
            print(f"  Min FDR: {summary_data.get('min_fdr', 1.0):.3f}")
            
            # Show top pathways
            top_pathways = summary_data.get('top_pathways', [])
            if top_pathways:
                print(f"  Top pathways:")
                for i, pathway in enumerate(top_pathways[:5]):  # Show top 5
                    term = pathway.get('Term', pathway.get('Name', 'Unknown'))
                    pval = pathway.get('P-value', pathway.get('FDR q-val', 1.0))
                    score = pathway.get('Enrichment Score', -np.log10(pval) if pval > 0 else 0)
                    print(f"    {i+1}. {term[:50]}... (score: {score:.3f}, p: {pval:.3e})")
        
        # Overall assessment
        print(f"\n" + "="*60)
        print("OVERALL ASSESSMENT")
        print("="*60)
        
        target_pathways_passed = summary.get('target_pathways_passed', 0)
        total_target_pathways = len(validation_results)
        
        if target_pathways_passed == total_target_pathways:
            print("✅ ALL TARGET PATHWAYS PASSED VALIDATION")
        elif target_pathways_passed > 0:
            print(f"⚠️  {target_pathways_passed}/{total_target_pathways} TARGET PATHWAYS PASSED")
        else:
            print("❌ NO TARGET PATHWAYS PASSED VALIDATION")
        
        print(f"\nAnalysis completed successfully!")
        print(f"Results saved to: pathway_impact_results.json")
        
        return True
        
    except Exception as e:
        logger.error(f"Error during pathway analysis: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main function"""
    print("Testing Pathway Impact Analysis with GSEApy")
    print("=" * 50)
    
    success = test_pathway_analysis()
    
    if success:
        print("\n✅ Test completed successfully!")
        return 0
    else:
        print("\n❌ Test failed!")
        return 1


if __name__ == "__main__":
    exit(main()) 