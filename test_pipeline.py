#!/usr/bin/env python3
"""
Test script for the Enhanced Pipeline Orchestrator

This script demonstrates how to use the enhanced pipeline orchestrator
with the example files provided in the project.
"""

import os
import sys
from pathlib import Path

# Add modules to path
sys.path.append(str(Path(__file__).parent / "modules"))

try:
    from enhanced_pipeline_orchestrator import EnhancedPipelineOrchestrator
except ImportError:
    print("Error: Could not import EnhancedPipelineOrchestrator")
    print("Please ensure the modules directory is properly set up")
    sys.exit(1)

def main():
    """Test the enhanced pipeline orchestrator with example files"""
    
    # Define input file paths
    variants_file = "examples/ex_variants.json"
    clinical_data_file = "examples/ex_clinical_data.json"
    protein_sequence_file = "examples/ex_protein_seq.txt"
    pdb_file = "examples/protein_42.pdb"
    
    # Check if all input files exist
    input_files = [variants_file, clinical_data_file, protein_sequence_file, pdb_file]
    missing_files = [f for f in input_files if not os.path.exists(f)]
    
    if missing_files:
        print(f"Error: Missing input files: {missing_files}")
        print("Please ensure all example files are available in the examples/ directory")
        return
    
    print("Enhanced Pipeline Orchestrator Test")
    print("=" * 40)
    print(f"Variants file: {variants_file}")
    print(f"Clinical data file: {clinical_data_file}")
    print(f"Protein sequence file: {protein_sequence_file}")
    print(f"PDB file: {pdb_file}")
    print()
    
    try:
        # Initialize orchestrator
        print("Initializing Enhanced Pipeline Orchestrator...")
        orchestrator = EnhancedPipelineOrchestrator(max_workers=4, output_dir="outputs")
        
        # Run pipeline
        print("Running pipeline...")
        result = orchestrator.run_pipeline(
            variants_file=variants_file,
            clinical_data_file=clinical_data_file,
            protein_sequence_file=protein_sequence_file,
            pdb_file=pdb_file
        )
        
        # Print results
        print("\nPipeline Results:")
        print("=" * 40)
        
        if result.success:
            print(f"✅ Pipeline completed successfully in {result.execution_time:.2f} seconds")
            print("\nOutput files generated:")
            print(f"  📄 Processed input: {result.processed_input_file}")
            print(f"  📄 Variant impact: {result.variant_impact_file}")
            print(f"  📄 Structural results: {result.structural_results_file}")
            print(f"  📄 Sequence variant: {result.sequence_variant_file}")
            print(f"  📄 Coherence results: {result.coherence_results_file}")
            print(f"  📄 Pathway dynamics: {result.pathway_dynamics_file}")
            
            # Check if output files exist
            output_files = [
                result.processed_input_file,
                result.variant_impact_file,
                result.structural_results_file,
                result.sequence_variant_file,
                result.coherence_results_file,
                result.pathway_dynamics_file
            ]
            
            existing_files = [f for f in output_files if os.path.exists(f)]
            print(f"\n📊 Generated {len(existing_files)}/{len(output_files)} output files")
            
        else:
            print(f"❌ Pipeline failed: {result.errors}")
            
        if result.warnings:
            print(f"\n⚠️  Warnings: {result.warnings}")
            
    except Exception as e:
        print(f"❌ Error running pipeline: {e}")
        import traceback
        traceback.print_exc()
        
    finally:
        # Cleanup
        try:
            orchestrator.cleanup()
        except Exception as e:
            print(f"Warning: Error during cleanup: {e}")


if __name__ == "__main__":
    main() 