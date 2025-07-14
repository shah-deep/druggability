#!/usr/bin/env python3
"""
Check variant file for null scores
"""

import json

def check_variants():
    with open('variant_impact_results_large.json', 'r') as f:
        data = json.load(f)
    
    variants = data['missense_variants']
    print(f"Total variants: {len(variants)}")
    
    # Check for null scores
    null_scores = [v for v in variants if v.get('pathogenicity_score') is None]
    print(f"Variants with null scores: {len(null_scores)}")
    
    if null_scores:
        print("Sample null score variants:")
        for v in null_scores[:5]:
            print(f"  {v['id']}: {v['gene']}")
    
    # Check unique genes
    genes = set(v['gene'] for v in variants if v.get('gene'))
    print(f"Unique genes: {len(genes)}")
    print(f"Sample genes: {list(genes)[:10]}")

if __name__ == "__main__":
    check_variants() 