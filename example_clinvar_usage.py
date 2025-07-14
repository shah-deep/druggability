#!/usr/bin/env python3
"""
Example usage of ClinVar annotation functionality

This script demonstrates how to integrate ClinVar annotation into your variant analysis workflow.
"""

import json
import sys
from pathlib import Path

# Add modules to path
sys.path.append(str(Path(__file__).parent / "modules"))

from clinvar_annotator import ClinVarAnnotator


def annotate_variants_with_clinvar(variants_data):
    """
    Annotate variants with ClinVar information
    
    Args:
        variants_data: List of variant dictionaries
        
    Returns:
        Dictionary with variant annotations
    """
    
    # Initialize ClinVar annotator
    annotator = ClinVarAnnotator()
    
    # Annotate all variants
    annotations = annotator.annotate_variants(variants_data)
    
    return annotations


def main():
    """Example usage with the provided variant data"""
    
    # Your variant data (from the user's example)
    variants = [
        {
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
    ]
    
    print("Annotating variants with ClinVar...")
    
    # Get ClinVar annotations
    annotations = annotate_variants_with_clinvar(variants)
    
    # Process results
    for variant_id, annotation in annotations.items():
        print(f"\nVariant: {variant_id}")
        print(f"  Gene: {annotation.gene}")
        print(f"  Protein Change: {annotation.protein_change}")
        print(f"  Clinical Significance: {annotation.clinical_significance}")
        
        # Example: Check if variant is pathogenic
        if annotation.clinical_significance:
            if "pathogenic" in annotation.clinical_significance.lower():
                print(f"  ⚠️  This variant is classified as PATHOGENIC")
            elif "benign" in annotation.clinical_significance.lower():
                print(f"  ✅ This variant is classified as BENIGN")
            else:
                print(f"  ❓ This variant has uncertain significance")
        
        print(f"  Review Status: {annotation.review_status}")
        print(f"  Condition: {annotation.condition}")
        print(f"  Source: {annotation.source}")
        
        if annotation.warnings:
            print(f"  Warnings: {', '.join(annotation.warnings)}")
    
    # Save results
    output_file = "clinvar_annotations.json"
    output_data = {}
    
    for variant_id, annotation in annotations.items():
        output_data[variant_id] = {
            "gene": annotation.gene,
            "protein_change": annotation.protein_change,
            "clinical_significance": annotation.clinical_significance,
            "review_status": annotation.review_status,
            "condition": annotation.condition,
            "variation_id": annotation.variation_id,
            "last_evaluated": annotation.last_evaluated,
            "submitter": annotation.submitter,
            "star_rating": annotation.star_rating,
            "accession": annotation.accession,
            "source": annotation.source,
            "warnings": annotation.warnings
        }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")
    
    return annotations


if __name__ == "__main__":
    main() 