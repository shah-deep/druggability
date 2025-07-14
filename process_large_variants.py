#!/usr/bin/env python3
"""
Script to process large_variants.json using enhanced_input_processor
"""

import json
import sys
import os
from pathlib import Path

# Add the modules directory to the path
sys.path.append(str(Path(__file__).parent / "modules"))

from enhanced_input_processor import EnhancedInputProcessor

def main():
    """Process the large_variants.json file using enhanced_input_processor"""
    
    # Load the variants data
    variants_file = "examples/large_variants.json"
    
    if not os.path.exists(variants_file):
        print(f"Error: {variants_file} not found")
        return
    
    print(f"Loading variants from {variants_file}...")
    
    with open(variants_file, 'r') as f:
        variants_data = json.load(f)
    
    print(f"Loaded {len(variants_data)} variants")
    
    # Create input dictionary for the processor
    inputs = {
        'missense_variants': variants_data,
        # Add a dummy PDB file to satisfy the required field
        'pdb_file': 'examples/protein_42.pdb' if os.path.exists('examples/protein_42.pdb') else None
    }
    
    # Initialize the enhanced input processor
    print("Initializing EnhancedInputProcessor...")
    processor = EnhancedInputProcessor()
    
    # Process the inputs
    print("Processing inputs...")
    processed_input = processor.process_inputs(inputs)
    
    # Get and display summary
    print("\n" + "="*50)
    print("PROCESSING SUMMARY")
    print("="*50)
    
    summary = processor.get_input_summary(processed_input)
    print(json.dumps(summary, indent=2))
    
    # Display validation results
    print("\n" + "="*50)
    print("VALIDATION RESULTS")
    print("="*50)
    
    if processed_input.validation_result:
        if processed_input.validation_result.is_valid:
            print("✅ Input validation PASSED")
        else:
            print("❌ Input validation FAILED")
            print("Errors:")
            for error in processed_input.validation_result.errors:
                print(f"  - {error}")
        
        if processed_input.validation_result.warnings:
            print("Warnings:")
            for warning in processed_input.validation_result.warnings:
                print(f"  - {warning}")
    else:
        print("⚠️  No validation result available")
    
    # Display variant statistics
    if processed_input.missense_variants:
        print("\n" + "="*50)
        print("VARIANT STATISTICS")
        print("="*50)
        
        total_variants = len(processed_input.missense_variants)
        print(f"Total variants: {total_variants}")
        
        # Count unique genes
        unique_genes = set()
        for variant in processed_input.missense_variants:
            if variant.get('gene'):
                unique_genes.add(variant['gene'])
        
        print(f"Unique genes: {len(unique_genes)}")
        
        # Show first few variants as example
        print("\nFirst 5 variants:")
        for i, variant in enumerate(processed_input.missense_variants[:5]):
            print(f"  {i+1}. {variant.get('gene', 'Unknown')} - {variant.get('protein_change', 'Unknown')}")
    
    # Export processed input if validation passed
    if processed_input.validation_result and processed_input.validation_result.is_valid:
        output_file = "processed_large_input.json"
        print(f"\nExporting processed input to {output_file}...")
        processor.export_processed_input(processed_input, output_file)
        print(f"✅ Successfully exported to {output_file}")
    else:
        print("\n❌ Skipping export due to validation failures")

if __name__ == "__main__":
    main() 