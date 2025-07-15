#!/usr/bin/env python3
"""
Enhanced Pipeline Orchestrator

Manages parallel processing workflows for variant analysis, structural analysis,
and pathway analysis. Coordinates multiple modules to process variants, clinical data,
protein sequences, and PDB structures.

Example usage:

To run the enhanced pipeline orchestrator, use the following command in WSL with your virtual environment activated:

    wsl bash -c "source .venv/bin/activate && python3 -m orchestrator.enhanced_pipeline_orchestrator \
        --variants examples/ex_variants.json \
        --clinical examples/ex_clinical_data.json \
        --protein-seq examples/ex_protein_seq.txt \
        --pdb examples/protein_42.pdb \
        --output-dir outputs \
        --max-workers 8"

"""

import os
import json
import logging
import multiprocessing as mp
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import tempfile
import shutil
import concurrent.futures
import asyncio
import argparse
from modules.enhanced_input_processor import EnhancedInputProcessor
from modules.variant_impact_analyzer import VariantImpactAnalyzer
from modules.structure_function_integrator import StructureFunctionIntegrator
from modules.sequence_variant_applier import SequenceVariantApplier
from modules.coherence_analyzer import CoherenceAnalyzer
from modules.pathway_dynamics_analyzer import PathwayImpactAnalyzer
from modules.intelligent_coherence_scorer import IntelligentCoherenceScorer

# Configure logging
from modules.logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)

# Set multiprocessing start method to 'spawn' to avoid daemon process issues
# Only set if not already set
if not hasattr(mp, '_start_method'):
    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        # If spawn is not available, use default
        pass


@dataclass
class PipelineResult:
    """Result of the complete pipeline execution"""
    processed_input_file: str
    variant_impact_file: str
    structural_results_file: str
    sequence_variant_file: str
    coherence_results_file: str
    pathway_dynamics_file: str
    final_results_file: str
    execution_time: float
    success: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class EnhancedPipelineOrchestrator:
    """
    Orchestrates parallel processing workflows for comprehensive variant analysis
    """
    
    def __init__(self, max_workers: Optional[int] = None, output_dir: str = "outputs"):
        """
        Initialize the pipeline orchestrator
        
        Args:
            max_workers: Maximum number of parallel processes (default: CPU count)
            output_dir: Output directory for results (default: "outputs")
        """
        self.max_workers = max_workers or mp.cpu_count()
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize modules
        self.input_processor = EnhancedInputProcessor()
        self.variant_analyzer = VariantImpactAnalyzer()
        self.structure_integrator = StructureFunctionIntegrator()
        self.sequence_applier = SequenceVariantApplier()
        self.coherence_analyzer = CoherenceAnalyzer()
        self.intelligent_scorer = IntelligentCoherenceScorer()
        
        
        # Mark pipeline start in all log files
        logger.info(f"=== Enhanced Pipeline Orchestrator START === {self.timestamp}")
        logger.warning(f"=== Enhanced Pipeline Orchestrator START === {self.timestamp}")
        logger.error(f"=== Enhanced Pipeline Orchestrator START === {self.timestamp}")
        logger.info(f"Orchestrator initialized with {self.max_workers} workers")

    
    def run_pipeline(self, 
                    variants_file: str,
                    clinical_data_file: str,
                    protein_sequence_file: str,
                    pdb_file: str) -> PipelineResult:
        """
        Run the complete pipeline with parallel processing
        
        Args:
            variants_file: Path to variants JSON file
            clinical_data_file: Path to clinical data JSON file
            protein_sequence_file: Path to protein sequence file
            pdb_file: Path to PDB structure file
            
        Returns:
            PipelineResult with all output file paths and execution status
        """
        start_time = time.time()
        result = PipelineResult(
            processed_input_file="",
            variant_impact_file="",
            structural_results_file="",
            sequence_variant_file="",
            coherence_results_file="",
            pathway_dynamics_file="",
            final_results_file="",
            execution_time=0.0,
            success=False
        )
        
        try:
            logger.info("Starting enhanced pipeline execution")
            
            # Step 1 & 2: Run initial analysis (input processing + variant impact) and structural analysis in parallel
            logger.info("Running initial analysis and structural analysis in parallel")
            
            # Prepare tasks for parallel execution
            with concurrent.futures.ProcessPoolExecutor(max_workers=min(self.max_workers, 2)) as executor:
                # Submit both tasks
                initial_future = executor.submit(self._run_initial_analysis_parallel, 
                                              variants_file, clinical_data_file, protein_sequence_file, pdb_file)
                structural_future = executor.submit(self._run_structural_analysis, pdb_file)
                
                # Wait for both to complete
                try:
                    initial_result = initial_future.result()
                    structural_result = structural_future.result()
                except Exception as e:
                    logger.error(f"Error in parallel execution: {e}")
                    result.errors.append(f"Parallel execution failed: {e}")
                    result.execution_time = time.time() - start_time
                    return result
            
            # Check for errors in parallel results
            if initial_result.get('type') == 'error':
                error_msg = initial_result.get('error', 'Unknown error in initial analysis')
                logger.error(f"Initial analysis failed: {error_msg}")
                result.errors.append(f"Initial analysis failed: {error_msg}")
                result.execution_time = time.time() - start_time
                return result
                
            if structural_result.get('type') == 'error':
                error_msg = structural_result.get('error', 'Unknown error in structural analysis')
                logger.error(f"Structural analysis failed: {error_msg}")
                result.errors.append(f"Structural analysis failed: {error_msg}")
                result.execution_time = time.time() - start_time
                return result
            
            # Collect results from parallel initial steps
            if initial_result['type'] == 'initial_analysis':
                result.processed_input_file = initial_result['processed_input_file']
                result.variant_impact_file = initial_result['variant_impact_file']
            elif initial_result['type'] == 'structural_analysis':
                result.structural_results_file = initial_result['file']
                
            if structural_result['type'] == 'structural_analysis':
                result.structural_results_file = structural_result['file']
            
            logger.info("Initial parallel steps complete. Starting secondary parallel processing...")
            
            # Steps 3-5: Run parallel processes
            try:
                parallel_results = self._run_parallel_processes(
                    result.processed_input_file, result.variant_impact_file, clinical_data_file, protein_sequence_file
                )
                
                result.sequence_variant_file = parallel_results['sequence_variant_file']
                result.coherence_results_file = parallel_results['coherence_results_file']
                result.pathway_dynamics_file = parallel_results['pathway_dynamics_file']
            except Exception as e:
                logger.error(f"Error in parallel processes: {e}")
                result.errors.append(f"Parallel processes failed: {e}")
                result.execution_time = time.time() - start_time
                return result
            
            # Step 6: Run intelligent coherence scoring
            try:
                final_results_file = self._run_intelligent_coherence_scoring(
                    result.structural_results_file,
                    result.variant_impact_file,
                    result.sequence_variant_file,
                    result.coherence_results_file,
                    result.pathway_dynamics_file
                )
                result.final_results_file = final_results_file
            except Exception as e:
                logger.error(f"Error in intelligent coherence scoring: {e}")
                result.errors.append(f"Intelligent coherence scoring failed: {e}")
                result.execution_time = time.time() - start_time
                return result
            
            result.success = True
            result.execution_time = time.time() - start_time
            
            logger.info(f"Pipeline execution completed successfully in {result.execution_time:.2f} seconds")
            
        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            result.execution_time = time.time() - start_time
            logger.error(f"Pipeline execution failed: {e}")
        
        return result
    
    
    def _process_inputs(self, 
                       variants_file: str,
                       clinical_data_file: str,
                       protein_sequence_file: str,
                       pdb_file: str) -> str:
        """Process all inputs using enhanced_input_processor"""
        logger.info("Processing inputs with EnhancedInputProcessor")
        
        # Load input data
        inputs = {}
        
        # Load variants
        if os.path.exists(variants_file):
            with open(variants_file, 'r') as f:
                inputs['missense_variants'] = json.load(f)
        
        # Load clinical data
        if os.path.exists(clinical_data_file):
            with open(clinical_data_file, 'r') as f:
                inputs['clinical_data'] = json.load(f)
        
        # # Load protein sequence
        # if os.path.exists(protein_sequence_file):
        #     with open(protein_sequence_file, 'r') as f:
        #         inputs['protein_sequence'] = f.read().strip()
        
        # Add PDB file path
        if os.path.exists(pdb_file):
            inputs['pdb_file'] = pdb_file
        
        # Process inputs
        processed_input = self.input_processor.process_inputs(inputs)
        
        # Save processed input
        output_file = self.output_dir / f"processed_input_{self.timestamp}.json"
        self.input_processor.export_processed_input(processed_input, str(output_file))
        
        logger.info(f"Processed inputs saved to: {output_file}")
        return str(output_file)
    
    def _run_variant_impact_analysis(self, processed_input_file: str) -> str:
        """Run variant impact analysis"""
        logger.info("Running variant impact analysis")
        
        output_file = self.output_dir / f"variant_impact_results_{self.timestamp}.json"
        
        # Run variant impact analysis
        results = self.variant_analyzer.analyze_variants(processed_input_file)
        self.variant_analyzer.save_results(results, str(output_file))
        
        logger.info(f"Variant impact analysis completed: {output_file}")
        return str(output_file)
    
    def _run_structural_analysis(self, pdb_file: str) -> Dict[str, str]:
        """Run structural analysis using structure_function_integrator"""
        logger.info(f"Running structural analysis for PDB file: {pdb_file}")
        
        # Validate PDB file exists and is readable
        if not os.path.exists(pdb_file):
            error_msg = f"PDB file not found: {pdb_file}"
            logger.error(error_msg)
            return {'type': 'error', 'file': '', 'error': error_msg}
        
        try:
            # Check file size
            file_size = os.path.getsize(pdb_file)
            logger.info(f"PDB file size: {file_size} bytes")
            
            if file_size == 0:
                error_msg = f"PDB file is empty: {pdb_file}"
                logger.error(error_msg)
                return {'type': 'error', 'file': '', 'error': error_msg}
            
            output_file = self.output_dir / f"structural_results_{self.timestamp}.json"
            
            # Run structural analysis
            result = self.structure_integrator.analyze_structure(pdb_file)
            self.structure_integrator.export_results(result, str(output_file))
            
            logger.info(f"Structural analysis completed: {output_file}")
            return {'type': 'structural_analysis', 'file': str(output_file)}
            
        except Exception as e:
            error_msg = f"Structural analysis failed: {e}"
            logger.error(error_msg)
            return {'type': 'error', 'file': '', 'error': error_msg}
    
    def _run_parallel_processes(self, 
                               processed_input_file: str,
                               variant_impact_file: str,
                               clinical_data_file: str,
                               protein_sequence_file: str) -> Dict[str, str]:
        """
        Run parallel processes for sequence variant application, coherence analysis, and pathway dynamics
        
        Returns:
            Dictionary with output file paths
        """
        logger.info("Starting parallel processes")
    
        # Prepare tasks for parallel execution
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(self.max_workers, 3)) as executor:
            # Submit all three tasks
            sequence_future = executor.submit(self._run_sequence_variant_analysis, processed_input_file, protein_sequence_file)
            coherence_future = executor.submit(self._run_coherence_analysis, variant_impact_file, clinical_data_file)
            pathway_future = executor.submit(self._run_pathway_dynamics_analysis, variant_impact_file)
            
            # Wait for all to complete and collect results
            results = {}
            futures = {
                'sequence_variant': sequence_future,
                'coherence': coherence_future,
                'pathway_dynamics': pathway_future
            }
            
            for name, future in futures.items():
                try:
                    result = future.result()
                    if result.get('type') == 'error':
                        logger.error(f"{name} failed: {result.get('error', 'Unknown error')}")
                        raise Exception(f"{name} failed: {result.get('error', 'Unknown error')}")
                    results[name] = result
                except Exception as e:
                    logger.error(f"Error in {name} process: {e}")
                    raise Exception(f"{name} process failed: {e}")
        
        # Collect results
        output_files = {}
        if results.get('sequence_variant', {}).get('type') == 'sequence_variant':
            output_files['sequence_variant_file'] = results['sequence_variant']['file']
        if results.get('coherence', {}).get('type') == 'coherence':
            output_files['coherence_results_file'] = results['coherence']['file']
        if results.get('pathway_dynamics', {}).get('type') == 'pathway_dynamics':
            output_files['pathway_dynamics_file'] = results['pathway_dynamics']['file']
        
        logger.info("Parallel processes completed")
        return output_files
    

    
    def _run_initial_analysis_parallel(self, 
                                     variants_file: str,
                                     clinical_data_file: str,
                                     protein_sequence_file: str,
                                     pdb_file: str) -> Dict[str, str]:
        """Run initial analysis in parallel processing context"""
        try:
            # Step 1: Process inputs
            processed_input_file = self._process_inputs(variants_file, clinical_data_file, protein_sequence_file, pdb_file)
            
            # Step 2: Run variant impact analysis
            variant_impact_file = self._run_variant_impact_analysis(processed_input_file)
            
            return {
                'type': 'initial_analysis',
                'processed_input_file': processed_input_file,
                'variant_impact_file': variant_impact_file
            }
        except Exception as e:
            logger.error(f"Error in initial analysis parallel process: {e}")
            return {'type': 'error', 'file': '', 'error': str(e)}
    
    def _run_sequence_variant_analysis(self, processed_input_file: str, protein_sequence_file: str) -> Dict[str, str]:
        """Run sequence variant analysis"""
        logger.info("Running sequence variant analysis")
        
        output_file = self.output_dir / f"sequence_variant_results_{self.timestamp}.json"
        
        # Load protein sequence if provided
        if os.path.exists(protein_sequence_file):
            with open(processed_input_file, 'r') as f:
                input_data = json.load(f)
            
            with open(protein_sequence_file, 'r') as f:
                protein_sequence = f.read().strip()
            
            input_data['protein_sequence'] = protein_sequence
            
            # Save updated input data
            temp_file = self.output_dir / f"temp_input_{self.timestamp}.json"
            with open(temp_file, 'w') as f:
                json.dump(input_data, f, indent=2)
            
            # Run sequence variant analysis
            results = self.sequence_applier.apply_variants_to_sequences(str(temp_file))
            self.sequence_applier.save_results(results, str(output_file))
            
            # Clean up temp file
            temp_file.unlink(missing_ok=True)
        else:
            # Run without protein sequence
            results = self.sequence_applier.apply_variants_to_sequences(processed_input_file)
            self.sequence_applier.save_results(results, str(output_file))
        
        logger.info(f"Sequence variant analysis completed: {output_file}")
        return {'type': 'sequence_variant', 'file': str(output_file)}
    
    def _run_coherence_analysis(self, variant_impact_file: str, clinical_data_file: str) -> Dict[str, str]:
        """Run coherence analysis"""
        logger.info("Running coherence analysis")
        
        output_file = self.output_dir / f"coherence_results_{self.timestamp}.json"
        
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        try:
            result = loop.run_until_complete(
                self.coherence_analyzer.analyze_coherence(variant_impact_file, clinical_data_file)
            )
            self.coherence_analyzer.export_results(result, str(output_file))
        finally:
            loop.close()
        
        logger.info(f"Coherence analysis completed: {output_file}")
        return {'type': 'coherence', 'file': str(output_file)}
    
    def _run_pathway_dynamics_analysis(self, variant_impact_file: str) -> Dict[str, str]:
        """Run pathway dynamics analysis"""
        logger.info("Running pathway dynamics analysis")
        
        output_file = self.output_dir / f"pathway_dynamics_results_{self.timestamp}.json"
        
        # Run pathway dynamics analysis
        analyzer = PathwayImpactAnalyzer(variant_impact_file)
        results = analyzer.run_complete_analysis(str(output_file))
        
        logger.info(f"Pathway dynamics analysis completed: {output_file}")
        return {'type': 'pathway_dynamics', 'file': str(output_file)}
    
    def _run_intelligent_coherence_scoring(self,
                                         structural_results_file: str,
                                         variant_impact_file: str,
                                         sequence_variant_file: str,
                                         coherence_results_file: str,
                                         pathway_dynamics_file: str) -> str:
        """Run intelligent coherence scoring"""
        logger.info("Running intelligent coherence scoring")
        
        output_file = self.output_dir / f"final_results_{self.timestamp}.json"
        
        # Run intelligent coherence scoring
        self.intelligent_scorer.run_scoring_pipeline(
            structural_results_file,
            variant_impact_file,
            sequence_variant_file,
            coherence_results_file,
            pathway_dynamics_file,
            str(output_file)
        )
        
        logger.info(f"Intelligent coherence scoring completed: {output_file}")
        return str(output_file)
    
    def cleanup(self):
        """Clean up temporary files and resources"""
        try:
            self.structure_integrator.cleanup()
        except Exception as e:
            logger.error(f"Error during cleanup: {e}")


def main():
    """Main function for command-line usage"""
    
    parser = argparse.ArgumentParser(description="Enhanced Pipeline Orchestrator")
    parser.add_argument("--variants", required=True, help="Path to variants JSON file")
    parser.add_argument("--clinical", required=True, help="Path to clinical data JSON file")
    parser.add_argument("--protein-seq", required=True, help="Path to protein sequence file")
    parser.add_argument("--pdb", required=True, help="Path to PDB structure file")
    parser.add_argument("--output-dir", default="outputs", help="Output directory for results (default: outputs)")
    parser.add_argument("--max-workers", type=int, help="Number of parallel workers (default: CPU count)")
    
    args = parser.parse_args()
    
    # Initialize orchestrator
    orchestrator = EnhancedPipelineOrchestrator(
        max_workers=args.max_workers,
        output_dir=args.output_dir
    )
    
    try:
        # Run pipeline
        result = orchestrator.run_pipeline(
            variants_file=args.variants,
            clinical_data_file=args.clinical,
            protein_sequence_file=args.protein_seq,
            pdb_file=args.pdb
        )
        
        # Print results
        if result.success:
            logger.info(f"Pipeline completed successfully in {result.execution_time:.2f} seconds (identifier: {result.timestamp})")
            logger.info(f"Output files:")
            logger.info(f"  Processed input: {result.processed_input_file}")
            logger.info(f"  Variant impact: {result.variant_impact_file}")
            logger.info(f"  Structural results: {result.structural_results_file}")
            logger.info(f"  Sequence variant: {result.sequence_variant_file}")
            logger.info(f"  Coherence results: {result.coherence_results_file}")
            logger.info(f"  Pathway dynamics: {result.pathway_dynamics_file}")
            logger.info(f"  Final results: {result.final_results_file}")
        else:
            logger.error(f"Pipeline failed: {result.errors} (identifier: {result.timestamp})")
            
    finally:
        orchestrator.cleanup()


if __name__ == "__main__":
    main() 