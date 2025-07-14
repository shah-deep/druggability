#!/usr/bin/env python3
"""
Test script for ClinVar annotation functionality
"""

import json
import sys
from pathlib import Path

# Add modules to path
sys.path.append(str(Path(__file__).parent / "modules"))

from clinvar_annotator import ClinVarAnnotator


def test_clinvar_annotation():
    """Test ClinVar annotation with the provided example variant"""
    
    # Example variant from the user
    test_variant = {
        "id": "var_001",
        "position": 7675088,
        "reference": "G",
        "alternate": "A",
        "gene": "TP53",
        "protein_change": "p.Arg175His",
        "chromosome": "17",
        "transcript": "",
        "transcript_id": "ENST00000269305",
        "vep_transcript_ids": [
            "ENST00000269305",
            "ENST00000359597",
            "ENST00000413465",
            "ENST00000420246",
            "ENST00000445888",
            "ENST00000455263",
            "ENST00000503591",
            "ENST00000504290",
            "ENST00000504937",
            "ENST00000505014",
            "ENST00000508793",
            "ENST00000509690",
            "ENST00000510385",
            "ENST00000514944",
            "ENST00000574684",
            "ENST00000576024",
            "ENST00000604348",
            "ENST00000610292",
            "ENST00000610538",
            "ENST00000610623",
            "ENST00000618944",
            "ENST00000619186",
            "ENST00000619485",
            "ENST00000620739",
            "ENST00000622645",
            "ENST00000635293",
            "ENST00000714356",
            "ENST00000714357",
            "ENST00000714358",
            "ENST00000714359",
            "ENST00000714408",
            "ENST00000714409"
        ]
    }
    
    # Additional test variants
    additional_variants = [
        {
            "id": "var_002",
            "position": 7674220,
            "reference": "T",
            "alternate": "C",
            "gene": "TP53",
            "protein_change": "p.Arg273His",
            "chromosome": "17",
            "transcript_id": "ENST00000269305"
        },
        {
            "id": "var_003",
            "position": 117199644,
            "reference": "T",
            "alternate": "C",
            "gene": "CFTR",
            "protein_change": "p.Phe508del",
            "chromosome": "7",
            "transcript_id": "ENST00000003084"
        }
    ]
    
    # Initialize annotator
    print("Initializing ClinVar Annotator...")
    annotator = ClinVarAnnotator()
    
    # Test single variant annotation
    print(f"\n{'='*80}")
    print("TESTING SINGLE VARIANT ANNOTATION")
    print(f"{'='*80}")
    
    result = annotator.annotate_variant(test_variant)
    
    print(f"Variant ID: {result.variant_id}")
    print(f"Gene: {result.gene}")
    print(f"Protein Change: {result.protein_change}")
    print(f"Clinical Significance: {result.clinical_significance}")
    print(f"Review Status: {result.review_status}")
    print(f"Condition: {result.condition}")
    print(f"Variation ID: {result.variation_id}")
    print(f"Last Evaluated: {result.last_evaluated}")
    print(f"Submitter: {result.submitter}")
    print(f"Star Rating: {result.star_rating}")
    print(f"Accession: {result.accession}")
    print(f"Source: {result.source}")
    
    if result.warnings:
        print(f"\nWarnings:")
        for warning in result.warnings:
            print(f"  - {warning}")
    
    # Test multiple variant annotation
    print(f"\n{'='*80}")
    print("TESTING MULTIPLE VARIANT ANNOTATION")
    print(f"{'='*80}")
    
    all_variants = [test_variant] + additional_variants
    results = annotator.annotate_variants(all_variants)
    
    for variant_id, result in results.items():
        print(f"\nVariant: {variant_id}")
        print(f"  Gene: {result.gene}")
        print(f"  Protein Change: {result.protein_change}")
        print(f"  Clinical Significance: {result.clinical_significance}")
        print(f"  Source: {result.source}")
        if result.warnings:
            print(f"  Warnings: {', '.join(result.warnings)}")
    
    # Save results to JSON
    output_file = "clinvar_annotation_results.json"
    output_data = {}
    
    for variant_id, result in results.items():
        output_data[variant_id] = {
            "gene": result.gene,
            "protein_change": result.protein_change,
            "clinical_significance": result.clinical_significance,
            "review_status": result.review_status,
            "condition": result.condition,
            "variation_id": result.variation_id,
            "last_evaluated": result.last_evaluated,
            "submitter": result.submitter,
            "star_rating": result.star_rating,
            "accession": result.accession,
            "source": result.source,
            "warnings": result.warnings
        }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")
    print(f"\n{'='*80}")


if __name__ == "__main__":
    test_clinvar_annotation() 