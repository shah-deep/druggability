#!/usr/bin/env python3
"""
Test script for Pathway Impact Analysis with GSEApy

This script runs the pathway impact analyzer on the processed_input.json file
and validates the minimum expected scores for target pathways.
"""

import sys
import os
import json
from datetime import datetime
from modules.variant_impact_analyzer import VariantImpactAnalyzer

# Add modules directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'modules'))


# Example variant from user
example_variant = {
    "id": "var_BRCA1_185delAG",
    "position": 43124,
    "reference": "AG",
    "alternate": "A",
    "gene": "BRCA1",
    "protein_change": "p.Glu23fs",
    "chromosome": "17",
    "transcript": "",
    "transcript_id": "ENST00000357654",
    "vep_transcript_ids": [
        "ENST00000497488", "ENST00000489037", "ENST00000478531", "ENST00000357654", "ENST00000473961",
        "ENST00000477152", "ENST00000352993", "ENST00000493919", "ENST00000494123", "ENST00000471181",
        "ENST00000652672", "ENST00000634433", "ENST00000476777", "ENST00000700081", "ENST00000470026",
        "ENST00000713676", "ENST00000618469", "ENST00000461574", "ENST00000644555", "ENST00000468300",
        "ENST00000700082", "ENST00000644379", "ENST00000484087", "ENST00000586385", "ENST00000591534",
        "ENST00000591849", "ENST00000493795", "ENST00000461221", "ENST00000491747", "ENST00000472490",
        "ENST00000700182", "ENST00000621897", "ENST00000354071", "ENST00000700183", "ENST00000492859",
        "ENST00000642945", "ENST00000700184", "ENST00000461798", "ENST00000700083", "ENST00000700185",
        "ENST00000700186"
    ]
}

def test_alphamissense_filtering():
    analyzer = VariantImpactAnalyzer()
    result = analyzer._analyze_alphamissense(
        example_variant["id"],
        example_variant["gene"],
        example_variant["protein_change"],
        example_variant["vep_transcript_ids"]
    )
    print("Total matches:", result.total_matches)
    print("Matching transcripts:", result.matching_transcripts)
    print("Warnings:", result.warnings)
    print("Average pathogenicity score:", result.average_pathogenicity_score)
    print("Max occurring class:", result.max_occurring_class)
    # Optionally, add asserts for expected behavior


def main():
    test_alphamissense_filtering()


if __name__ == "__main__":
    main() 