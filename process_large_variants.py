#!/usr/bin/env python3
"""
Process large variants through enhanced input processor and variant impact analyzer
"""

import json
import os
import sys
from pathlib import Path

# Add modules to path
sys.path.append('modules')

from enhanced_input_processor import EnhancedInputProcessor
from variant_impact_analyzer import VariantImpactAnalyzer

def main():
    """Process large variants and generate new variant impact results"""
    
    print("Loading large variant file...")
    with open('large_variants.json', 'r') as f:
        variants = json.load(f)
    
    print(f"Loaded {len(variants)} variants")
    
    # Create input for enhanced processor
    input_data = {
        'pdb_file': 'examples/protein_42.pdb',
        'missense_variants': variants
    }
    
    print("Processing through enhanced input processor...")
    processor = EnhancedInputProcessor()
    processed_input = processor.process_inputs(input_data)
    
    # Save processed input
    with open('processed_large_input.json', 'w') as f:
        json.dump(processed_input.__dict__, f, indent=2, default=str)
    
    print("Processed input saved to processed_large_input.json")
    
    # Apply to variant impact analyzer
    print("Applying to variant impact analyzer...")
    analyzer = VariantImpactAnalyzer()
    results = analyzer.analyze_variants('processed_large_input.json')
    
    # Save results
    with open('large_variant_impact_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("Results saved to large_variant_impact_results.json")
    
    # Print summary
    variants_with_scores = [v for v in results['missense_variants'] if v.get('pathogenicity_score')]
    print(f"\nSummary:")
    print(f"Total variants processed: {len(results['missense_variants'])}")
    print(f"Variants with pathogenicity scores: {len(variants_with_scores)}")
    
    if variants_with_scores:
        scores = [v['pathogenicity_score'] for v in variants_with_scores]
        print(f"Average pathogenicity score: {sum(scores)/len(scores):.4f}")
        print(f"Min pathogenicity score: {min(scores):.4f}")
        print(f"Max pathogenicity score: {max(scores):.4f}")
    
    # Count genes
    genes = set(v['gene'] for v in results['missense_variants'])
    print(f"Unique genes: {len(genes)}")
    print(f"Genes: {', '.join(sorted(genes))}")

if __name__ == "__main__":
    main() 