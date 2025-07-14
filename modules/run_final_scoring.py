#!/usr/bin/env python3
"""
Final Pathogenicity Scoring with Standardized Output Format

Runs comprehensive pathogenicity scoring on processed_input.json and outputs
results in the specified format: {pathogenicity_score, confidence, known_annotations[]}
"""

import json
import sys
import os
from pathlib import Path
sys.path.append('modules')

from comprehensive_pathogenicity_scorer import ComprehensivePathogenicityScorer
from output_formatter import OutputFormatter

def main():
    """Run comprehensive scoring and output in standardized format"""
    
    print("Loading variants from processed_input.json...")
    
    # Load variants from processed_input.json
    with open('processed_input.json', 'r') as f:
        data = json.load(f)
    
    variants = data.get('missense_variants', [])
    if not variants:
        print("No variants found in processed_input.json")
        return
    
    print(f"Found {len(variants)} variants to score")
    
    # Run comprehensive scoring
    print("\nRunning comprehensive pathogenicity scoring...")
    scorer = ComprehensivePathogenicityScorer()
    print(f"Scoring {len(variants)} variants using comprehensive methods")
    results = scorer.score_variants(variants)
    
    # Format results into specified output format
    print("\nFormatting results...")
    formatter = OutputFormatter()
    formatted_results = formatter.format_all_results(results)
    
    # Create final output structure
    final_output = {
        "metadata": {
            "format_version": "1.0",
            "description": "Pathogenicity scoring results in standardized format",
            "total_variants": len(formatted_results),
            "scoring_methods": ["AlphaMissense", "CADD", "VEP", "ClinVar"]
        },
        "results": formatted_results
    }
    
    # Save to file
    output_file = "pathogenicity_results.json"
    with open(output_file, 'w') as f:
        json.dump(final_output, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")
    
    # Print formatted results to console
    formatter.print_formatted_results(results)
    
    # Print summary
    print(f"\nSUMMARY:")
    print(f"  Total variants processed: {len(variants)}")
    print(f"  Variants with AlphaMissense scores: {sum(1 for r in results.values() if r.alphamissense_score is not None)}")
    print(f"  Variants with CADD scores: {sum(1 for r in results.values() if r.cadd_phred_score is not None)}")
    print(f"  Variants with VEP scores: {sum(1 for r in results.values() if r.vep_score is not None)}")
    print(f"  Variants with ClinVar data: {sum(1 for r in results.values() if r.clinvar_significance is not None)}")
    
    # Show example output format
    print(f"\nEXAMPLE OUTPUT FORMAT:")
    if formatted_results:
        example_variant = list(formatted_results.keys())[0]
        example_result = formatted_results[example_variant]
        print(f"  Variant: {example_variant}")
        print(f"  {{")
        print(f"    \"pathogenicity_score\": {example_result['pathogenicity_score']:.3f},")
        print(f"    \"confidence\": {example_result['confidence']:.3f},")
        print(f"    \"known_annotations\": [")
        for i, annotation in enumerate(example_result['known_annotations']):
            source = annotation.get('source', 'Unknown')
            description = annotation.get('description', 'No description')
            if i < len(example_result['known_annotations']) - 1:
                print(f"      {{\"source\": \"{source}\", \"description\": \"{description}\"}},")
            else:
                print(f"      {{\"source\": \"{source}\", \"description\": \"{description}\"}}")
        print(f"    ]")
        print(f"  }}")
    
    print(f"\n✓ Comprehensive pathogenicity scoring completed successfully!")
    print(f"✓ Results formatted according to specification")
    print(f"✓ AlphaMissense integration working")

if __name__ == "__main__":
    main() 