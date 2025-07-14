import pytest
import os
import json
from modules.variant_impact_analyzer import VariantImpactAnalyzer

def test_variant_impact_analyzer_expected_results(tmp_path):
    # Use the real processed_input.json
    input_file = 'processed_input.json'
    analyzer = VariantImpactAnalyzer()
    results = analyzer.analyze_variants(input_file)
    
    # Map expected values from the table, with upper limit set to 1.0
    expected = {
        ("TP53", "p.Arg175His"): {"score_range": (0.85, 1.0), "clinvar": ["Pathogenic", "Likely_Pathogenic"]},
        ("CFTR", "p.Phe508del"): {"score_range": (0.80, 1.0), "clinvar": ["Pathogenic"]},
        ("BRCA1", "p.Glu23fs"): {"score_range": (0.1, 0.2), "clinvar": ["Pathogenic"]},
    }
    
    found = 0
    for variant in results["missense_variants"]:
        key = (variant["gene"], variant["protein_change"])
        if key in expected:
            found += 1
            score = variant["pathogenicity_score"]
            lo, hi = expected[key]["score_range"]
            assert lo <= score <= hi, f"{key} pathogenicity_score {score} not in {lo}-{hi}"
            valid_clinvars = [v.lower() for v in expected[key]["clinvar"]]
            assert variant["clinvar_annotation"].lower() in valid_clinvars, f"{key} clinvar_annotation {variant['clinvar_annotation']} not in {expected[key]['clinvar']}"
    assert found == 3, f"Expected 3 variants, found {found}" 