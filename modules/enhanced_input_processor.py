#!/usr/bin/env python3
"""
Enhanced Input Processor for ITS4.2 Mechanistic Coherence Module
Handles multi-modal input validation, preprocessing, and routing for mechanistic coherence analysis
"""

import os
import json
import numpy as np
# import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any
from dataclasses import dataclass, field
import logging
import hashlib
from datetime import datetime
import warnings
import requests
import multiprocessing as mp
import tempfile
import pickle
import shutil
try:
    from .transcript_resolver import TranscriptResolver
except ImportError:
    from transcript_resolver import TranscriptResolver

# Configure logging
from .logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)

@dataclass
class InputValidationResult:
    """Result of input validation"""
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

@dataclass
class ProcessedInput:
    """Processed input data structure"""
    # Core structural data
    pdb_file: Optional[str] = None
    pocket_data: Optional[Dict] = None
    
    # Sequence data
    dna_sequence: Optional[str] = None
    rna_sequence: Optional[str] = None
    protein_sequence: Optional[str] = None
    
    # Variant and clinical data
    missense_variants: Optional[List[Dict]] = None
    clinical_data: Optional[Dict] = None
    
    # Additional data types
    structural_constraints: Optional[Dict] = None
    single_cell_data: Optional[str] = None  # Path to h5ad file
    metabolic_data: Optional[Dict] = None
    
    # Processing metadata
    input_hash: Optional[str] = None
    processing_timestamp: Optional[str] = None
    validation_result: Optional[InputValidationResult] = None

class EnhancedInputProcessor:
    """
    Enhanced input processor for ITS4.2 mechanistic coherence analysis
    Handles validation, preprocessing, and routing of multi-modal inputs
    """
    
    def __init__(self, config: Optional[Dict] = None):
        self.config = config or self._get_default_config()
        self.logger = logging.getLogger(__name__)
        self.validation_rules = self._initialize_validation_rules()
        self.transcript_resolver = TranscriptResolver()
        # Set number of parallel workers for variant processing
        self.max_workers = min(mp.cpu_count(), 8)  # Cap at 8 workers
        
    def _get_default_config(self) -> Dict:
        """Get default configuration"""
        return {
            'max_sequence_length': 1000000,
            'allowed_file_types': ['.pdb', '.json', '.h5ad', '.csv', '.txt'],
            'required_fields': [],  # No mandatory fields - all are optional
            'optional_fields': [
                'pdb_file', 'dna_sequence', 'rna_sequence', 'protein_sequence',
                'missense_variants', 'clinical_data', 'single_cell_data'
            ],
            'validation_strict': True,
            # Model-specific requirements
            'model_requirements': {
                'genotype_models': {
                    'enformer': ['dna_sequence'],
                    'bigrna': ['rna_sequence'],
                    'alphamissense': ['missense_variants'],
                    'genegenie': ['dna_sequence']
                },
                'simulation_models': {
                    'scgpt': ['single_cell_data'],
                    'scfoundation': ['single_cell_data'],
                    'octo': ['clinical_data']
                },
                'structure_models': {
                    'proteinmpnn': ['protein_sequence'],
                    'rfdiffusion': ['pdb_file'],
                    'esm3': ['protein_sequence'],
                    'atom1': ['rna_sequence']
                },
                'biomarker_models': {
                    'arteraai': ['clinical_data']
                },
                'metabolic_models': {
                    'cobrapy': ['dna_sequence']
                }
            }
        }
    
    def _initialize_validation_rules(self) -> Dict:
        """Initialize validation rules for different input types"""
        return {
            'pdb_file': {
                'required': False,
                'file_exists': True,
                'file_extension': '.pdb',
                'max_size_mb': 100
            },
            'dna_sequence': {
                'required': False,
                'pattern': r'^[ATGCatgc\s]+\.{0,3}$',
                'min_length': 10,
                'max_length': 1000000,
                'allow_file_path': True
            },
            'rna_sequence': {
                'required': False,
                'pattern': r'^[AUGCaugc\s]+\.{0,3}$',
                'min_length': 10,
                'max_length': 1000000
            },
            'protein_sequence': {
                'required': False,
                'pattern': r'^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy\s]+$',
                'min_length': 10,
                'max_length': 100000
            },
            'missense_variants': {
                'required': False,
                'schema': {
                    'required_fields': ['id', 'position', 'reference', 'alternate'],
                    'optional_fields': ['gene', 'protein_change', 'chromosome', 'transcript']
                }
            },
            'clinical_data': {
                'required': False,
                'schema': {
                    'required_fields': ['patient_id'],
                    'optional_fields': ['age', 'sex', 'diagnosis', 'treatments']
                }
            }
        }
    
    def process_inputs(self, inputs: Dict[str, Any]) -> ProcessedInput:
        """
        Main processing function for inputs
        
        Args:
            inputs: Dictionary containing input data
            
        Returns:
            ProcessedInput object with validated and processed data
        """
        self.logger.info("Starting enhanced input processing")
        
        # Create processed input object
        processed_input = ProcessedInput()
        
        # Validate inputs
        validation_result = self._validate_inputs(inputs)
        processed_input.validation_result = validation_result
        
        if not validation_result.is_valid:
            self.logger.error(f"Input validation failed: {validation_result.errors}")
            return processed_input
        
        # Process and store inputs
        processed_input = self._process_core_inputs(inputs, processed_input)
        processed_input = self._process_sequence_inputs(inputs, processed_input)
        processed_input = self._process_variant_inputs(inputs, processed_input)
        processed_input = self._process_clinical_inputs(inputs, processed_input)
        processed_input = self._process_additional_inputs(inputs, processed_input)
        
        # Add metadata
        processed_input.input_hash = self._generate_input_hash(inputs)
        processed_input.processing_timestamp = datetime.now().isoformat()
        
        self.logger.info("Input processing completed successfully")
        return processed_input
    
    def _validate_inputs(self, inputs: Dict[str, Any]) -> InputValidationResult:
        """Validate all inputs according to rules"""
        result = InputValidationResult(is_valid=True)
        
        # Check if any input is provided
        if not inputs or all(value is None for value in inputs.values()):
            result.is_valid = False
            result.errors.append("At least one input field must be provided")
            return result
        
        # Check required fields
        for field in self.config['required_fields']:
            if field not in inputs or inputs[field] is None:
                result.is_valid = False
                result.errors.append(f"Required field '{field}' is missing")
        
        # Validate each input type
        for field, value in inputs.items():
            if value is None:
                continue
                
            if field in self.validation_rules:
                field_validation = self._validate_field(field, value)
                if not field_validation['is_valid']:
                    result.is_valid = False
                    result.errors.extend(field_validation['errors'])
                result.warnings.extend(field_validation['warnings'])
        
        # Validate model-specific requirements
        model_validation = self._validate_model_requirements(inputs)
        if not model_validation['is_valid']:
            result.is_valid = False
            result.errors.extend(model_validation['errors'])
        result.warnings.extend(model_validation['warnings'])
        
        return result
    
    def _validate_field(self, field: str, value: Any) -> Dict[str, Any]:
        """Validate a specific field"""
        rules = self.validation_rules.get(field, {})
        errors = []
        warnings = []
        
        # File validation
        if rules.get('file_exists', False):
            if not os.path.exists(str(value)):
                errors.append(f"File '{value}' does not exist")
        
        # Handle file paths for sequences
        if rules.get('allow_file_path', False) and isinstance(value, str):
            if os.path.exists(value):
                # It's a valid file path, skip pattern and length validation
                return {
                    'is_valid': True,
                    'errors': [],
                    'warnings': []
                }
        
        # Pattern validation for sequences
        if 'pattern' in rules:
            import re
            if not re.match(rules['pattern'], str(value)):
                errors.append(f"Field '{field}' does not match expected pattern")
        
        # Length validation
        if 'min_length' in rules and len(str(value)) < rules['min_length']:
            errors.append(f"Field '{field}' is too short (min: {rules['min_length']})")
        
        if 'max_length' in rules and len(str(value)) > rules['max_length']:
            errors.append(f"Field '{field}' is too long (max: {rules['max_length']})")
        
        # Schema validation for complex objects
        if 'schema' in rules and isinstance(value, (list, dict)):
            schema_errors = self._validate_schema(value, rules['schema'])
            errors.extend(schema_errors)
        
        return {
            'is_valid': len(errors) == 0,
            'errors': errors,
            'warnings': warnings
        }
    
    def _validate_schema(self, data: Any, schema: Dict) -> List[str]:
        """Validate data against schema"""
        errors = []
        
        if isinstance(data, list):
            for i, item in enumerate(data):
                if not isinstance(item, dict):
                    errors.append(f"Item {i} is not a dictionary")
                    continue
                    
                for required_field in schema.get('required_fields', []):
                    if required_field not in item:
                        errors.append(f"Item {i} missing required field '{required_field}'")
        
        elif isinstance(data, dict):
            for required_field in schema.get('required_fields', []):
                if required_field not in data:
                    errors.append(f"Missing required field '{required_field}'")
        
        return errors
    
    def _validate_model_requirements(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Validate model-specific requirements"""
        errors = []
        warnings = []
        
        # Check which models can be run with available inputs
        available_models = []
        missing_requirements = []
        
        model_requirements = self.config.get('model_requirements', {})
        
        for category, models in model_requirements.items():
            for model, required_fields in models.items():
                can_run = True
                missing_fields = []
                
                for field in required_fields:
                    if field not in inputs or inputs[field] is None:
                        can_run = False
                        missing_fields.append(field)
                
                if can_run:
                    available_models.append(model)
                else:
                    missing_requirements.append({
                        'model': model,
                        'category': category,
                        'missing_fields': missing_fields
                    })
        
        # Add warnings for missing model requirements
        for req in missing_requirements:
            warnings.append(
                f"Model '{req['model']}' ({req['category']}) cannot run: "
                f"missing {', '.join(req['missing_fields'])}"
            )
        
        # Add info about available models
        if available_models:
            self.logger.info(f"Available models: {', '.join(available_models)}")
        
        return {
            'is_valid': len(errors) == 0,
            'errors': errors,
            'warnings': warnings,
            'available_models': available_models,
            'missing_requirements': missing_requirements
        }
    
    def _process_core_inputs(self, inputs: Dict[str, Any], processed: ProcessedInput) -> ProcessedInput:
        """Process core structural inputs"""
        if 'pdb_file' in inputs:
            processed.pdb_file = str(inputs['pdb_file'])
            
        if 'pocket_data' in inputs:
            processed.pocket_data = inputs['pocket_data']
            
        return processed
    
    def _process_sequence_inputs(self, inputs: Dict[str, Any], processed: ProcessedInput) -> ProcessedInput:
        """Process sequence inputs"""
        if 'dna_sequence' in inputs:
            dna_input = inputs['dna_sequence']
            if isinstance(dna_input, str) and os.path.exists(dna_input):
                # If it's a file path, store the path instead of loading the entire sequence
                processed.dna_sequence = dna_input
            else:
                # If it's actual sequence data, normalize it
                processed.dna_sequence = self._normalize_sequence(dna_input)
            
        if 'rna_sequence' in inputs:
            rna_input = inputs['rna_sequence']
            if isinstance(rna_input, str) and os.path.exists(rna_input):
                # If it's a file path, store the path instead of loading the entire sequence
                processed.rna_sequence = rna_input
            else:
                # If it's actual sequence data, normalize it
                processed.rna_sequence = self._normalize_sequence(rna_input)
            
        if 'protein_sequence' in inputs:
            protein_input = inputs['protein_sequence']
            if isinstance(protein_input, str) and os.path.exists(protein_input):
                # If it's a file path, store the path instead of loading the entire sequence
                processed.protein_sequence = protein_input
            else:
                # If it's actual sequence data, normalize it
                processed.protein_sequence = self._normalize_sequence(protein_input)
            
        return processed
    
    def _process_variant_inputs(self, inputs: Dict[str, Any], processed: ProcessedInput) -> ProcessedInput:
        """Process variant inputs with parallel transcript resolution"""
        if 'missense_variants' in inputs:
            variants = inputs['missense_variants']
            if isinstance(variants, str):
                # Try to load from file
                try:
                    with open(variants, 'r') as f:
                        variants = json.load(f)
                except Exception as e:
                    self.logger.warning(f"Could not load variants from file: {e}")
                    variants = []
            
            normalized_variants = self._normalize_variants(variants)
            
            # Resolve transcript IDs for variants using parallel processing
            if normalized_variants:
                self.logger.info(f"Processing {len(normalized_variants)} variants with {self.max_workers} parallel workers")
                processed.missense_variants = self._resolve_transcript_ids_parallel(normalized_variants)
            else:
                processed.missense_variants = []
            
        return processed
    
    def _process_clinical_inputs(self, inputs: Dict[str, Any], processed: ProcessedInput) -> ProcessedInput:
        """Process clinical inputs"""
        if 'clinical_data' in inputs:
            clinical_data = inputs['clinical_data']
            if isinstance(clinical_data, str):
                # Try to load from file
                try:
                    with open(clinical_data, 'r') as f:
                        clinical_data = json.load(f)
                except Exception as e:
                    self.logger.warning(f"Could not load clinical data from file: {e}")
                    clinical_data = {}
            
            processed.clinical_data = clinical_data
            
        return processed
    
    def _process_additional_inputs(self, inputs: Dict[str, Any], processed: ProcessedInput) -> ProcessedInput:
        """Process additional input types"""
        if 'structural_constraints' in inputs:
            processed.structural_constraints = inputs['structural_constraints']
            
        if 'single_cell_data' in inputs:
            processed.single_cell_data = str(inputs['single_cell_data'])
            
        if 'metabolic_data' in inputs:
            processed.metabolic_data = inputs['metabolic_data']
            
        return processed
    
    def _normalize_sequence(self, sequence: str) -> str:
        """Normalize sequence data"""
        if not sequence:
            return ""
        
        # Remove whitespace and convert to uppercase
        normalized = sequence.strip().upper()
        
        # Remove any non-sequence characters
        import re
        if 'T' in normalized and 'U' in normalized:
            # Mixed DNA/RNA, keep as is
            pass
        elif 'U' in normalized:
            # RNA sequence
            normalized = re.sub(r'[^AUGC]', '', normalized)
        else:
            # DNA or protein sequence
            normalized = re.sub(r'[^ATGC]', '', normalized)
        
        return normalized
    
    def _normalize_variants(self, variants: List[Dict]) -> List[Dict]:
        """Normalize variant data"""
        normalized = []
        
        for variant in variants:
            if not isinstance(variant, dict):
                continue
                
            # Ensure required fields are present and preserve id field
            normalized_variant = {
                'id': str(variant.get('id', '')),
                'position': int(variant.get('position', 0)),
                'reference': str(variant.get('reference', '')),
                'alternate': str(variant.get('alternate', '')),
                'gene': str(variant.get('gene', '')),
                'protein_change': str(variant.get('protein_change', '')),
                'chromosome': str(variant.get('chromosome', '')),
                'transcript': str(variant.get('transcript', '')),
                'gencode_id': str(variant.get('gencode_id', ''))
            }
            
            normalized.append(normalized_variant)
        
        return normalized
    
    def _resolve_transcript_ids_parallel(self, variants: List[Dict]) -> List[Dict]:
        """
        Resolve transcript IDs for variants using parallel processing
        
        Args:
            variants: List of variant dictionaries
            
        Returns:
            List of variants with resolved transcript IDs
        """
        if not variants:
            return []
        
        # Create temporary directory for inter-process communication
        temp_dir = Path(tempfile.mkdtemp(prefix="variant_processing_"))
        
        try:
            # Prepare arguments for each process
            process_args = []
            for i, variant in enumerate(variants):
                # Create unique temp file for this variant
                temp_file = temp_dir / f"variant_{i}_{variant['id']}.pkl"
                process_args.append((variant, str(temp_file)))
            
            # Process variants in parallel
            with mp.Pool(processes=self.max_workers) as pool:
                pool.map(self._process_single_variant_worker, process_args)
            
            # Collect results from temp files
            resolved_variants = []
            for i, variant in enumerate(variants):
                temp_file = temp_dir / f"variant_{i}_{variant['id']}.pkl"
                if temp_file.exists():
                    try:
                        with open(temp_file, 'rb') as f:
                            resolved_variant = pickle.load(f)
                        resolved_variants.append(resolved_variant)
                    except Exception as e:
                        self.logger.error(f"Error loading result for {variant['id']}: {e}")
                        # Add original variant if loading fails
                        resolved_variants.append(variant)
                else:
                    # Add original variant if temp file doesn't exist
                    resolved_variants.append(variant)
            
            return resolved_variants
            
        finally:
            # Clean up temporary directory
            shutil.rmtree(temp_dir, ignore_errors=True)
            self.logger.info(f"Cleaned up temporary directory: {temp_dir}")

    @staticmethod
    def _process_single_variant_worker(args: Tuple[Dict, str]) -> None:
        """
        Worker function to process a single variant in a separate process
        
        Args:
            args: Tuple of (variant, temp_file_path)
        """
        variant, temp_file_path = args
        variant_id = variant['id']
        
        try:
            # Create transcript resolver instance for this process
            resolver = TranscriptResolver()
            
            # Resolve transcript ID for this variant using the correct method
            resolved_variant = EnhancedInputProcessor._resolve_single_variant_with_resolver(resolver, variant)
            
            # Save result to temp file
            with open(temp_file_path, 'wb') as f:
                pickle.dump(resolved_variant, f)
                
            logger.info(f"Processed variant {variant_id}")
            
        except Exception as e:
            logger.error(f"Error processing variant {variant_id}: {e}")
            # Save original variant if processing fails
            with open(temp_file_path, 'wb') as f:
                pickle.dump(variant, f)

    @staticmethod
    def _resolve_single_variant_with_resolver(resolver: TranscriptResolver, variant: Dict) -> Dict:
        """
        Resolve transcript ID for a single variant using the TranscriptResolver
        
        Args:
            resolver: TranscriptResolver instance
            variant: Variant dictionary
            
        Returns:
            Variant dictionary with resolved transcript ID and gencode_id
        """
        # Use the resolve_variants method which handles everything internally
        resolved_variants = resolver.resolve_variants([variant])
        
        # Get the first (and only) resolved variant
        if resolved_variants:
            resolved_variant = resolved_variants[0]
            
            # Add gencode_id if gene is available
            gene = resolved_variant.get('gene')
            if gene and not resolved_variant.get('gencode_id'):
                try:
                    gencode_id = resolver.get_gencode_id(gene)
                    if gencode_id:
                        resolved_variant['gencode_id'] = gencode_id
                        logger.info(f"Added gencode_id for {gene}: {gencode_id}")
                except Exception as e:
                    logger.warning(f"Could not fetch gencode_id for {gene}: {e}")
            
            return resolved_variant
        else:
            return variant
    
    def _resolve_transcript_ids(self, variants: List[Dict]) -> List[Dict]:
        """Resolve transcript IDs for variants using VEP API (sequential fallback)"""
        return self.transcript_resolver.resolve_variants(variants)
    
    def _generate_input_hash(self, inputs: Dict[str, Any]) -> str:
        """Generate hash for input data"""
        # Create a stable representation of inputs
        input_str = json.dumps(inputs, sort_keys=True, default=str)
        return hashlib.md5(input_str.encode()).hexdigest()
    
    def get_input_summary(self, processed_input: ProcessedInput) -> Dict[str, Any]:
        """Get summary of processed inputs"""
        summary = {
            'has_pdb': processed_input.pdb_file is not None,
            'has_pocket_data': processed_input.pocket_data is not None,
            'has_dna_sequence': processed_input.dna_sequence is not None,
            'has_rna_sequence': processed_input.rna_sequence is not None,
            'has_protein_sequence': processed_input.protein_sequence is not None,
            'has_variants': processed_input.missense_variants is not None,
            'has_clinical_data': processed_input.clinical_data is not None,
            'has_single_cell_data': processed_input.single_cell_data is not None,
            'input_hash': processed_input.input_hash,
            'processing_timestamp': processed_input.processing_timestamp,
            'validation_passed': processed_input.validation_result.is_valid if processed_input.validation_result else False,
            'parallel_workers': self.max_workers
        }
        
        # Add counts
        if processed_input.missense_variants:
            summary['variant_count'] = len(processed_input.missense_variants)
            
            # Add transcript resolution statistics
            variants_with_transcripts = sum(1 for v in processed_input.missense_variants if v.get('transcript_id'))
            summary['variants_with_transcript_ids'] = variants_with_transcripts
            summary['transcript_resolution_rate'] = variants_with_transcripts / len(processed_input.missense_variants) if processed_input.missense_variants else 0.0
            
            # Add gencode_id resolution statistics
            variants_with_gencode_ids = sum(1 for v in processed_input.missense_variants if v.get('gencode_id'))
            summary['variants_with_gencode_ids'] = variants_with_gencode_ids
            summary['gencode_id_resolution_rate'] = variants_with_gencode_ids / len(processed_input.missense_variants) if processed_input.missense_variants else 0.0
        
        if processed_input.dna_sequence:
            if os.path.exists(processed_input.dna_sequence):
                # It's a file path, get file size instead of sequence length
                try:
                    file_size = os.path.getsize(processed_input.dna_sequence)
                    summary['dna_file_size_bytes'] = file_size
                    summary['dna_file_path'] = processed_input.dna_sequence
                except OSError:
                    summary['dna_file_path'] = processed_input.dna_sequence
            else:
                # It's actual sequence data
                summary['dna_length'] = len(processed_input.dna_sequence)
        
        if processed_input.rna_sequence:
            summary['rna_length'] = len(processed_input.rna_sequence)
        
        if processed_input.protein_sequence:
            summary['protein_length'] = len(processed_input.protein_sequence)
        
        return summary
    
    def export_processed_input(self, processed_input: ProcessedInput, output_path: str) -> None:
        """Export processed input to file"""
        output_data = {}
        
        # Only include fields that have actual values
        if processed_input.pdb_file is not None:
            output_data['pdb_file'] = processed_input.pdb_file
        if processed_input.pocket_data is not None:
            output_data['pocket_data'] = processed_input.pocket_data
        if processed_input.dna_sequence is not None:
            output_data['dna_sequence'] = processed_input.dna_sequence
        if processed_input.rna_sequence is not None:
            output_data['rna_sequence'] = processed_input.rna_sequence
        if processed_input.protein_sequence is not None:
            output_data['protein_sequence'] = processed_input.protein_sequence
        if processed_input.missense_variants is not None:
            output_data['missense_variants'] = processed_input.missense_variants
        if processed_input.clinical_data is not None:
            output_data['clinical_data'] = processed_input.clinical_data
        if processed_input.structural_constraints is not None:
            output_data['structural_constraints'] = processed_input.structural_constraints
        if processed_input.single_cell_data is not None:
            output_data['single_cell_data'] = processed_input.single_cell_data
        if processed_input.metabolic_data is not None:
            output_data['metabolic_data'] = processed_input.metabolic_data
        
        # Always include metadata
        output_data['input_hash'] = processed_input.input_hash
        output_data['processing_timestamp'] = processed_input.processing_timestamp
        output_data['validation_result'] = {
            'is_valid': processed_input.validation_result.is_valid,
            'errors': processed_input.validation_result.errors,
            'warnings': processed_input.validation_result.warnings
        } if processed_input.validation_result else None
        output_data['parallel_processing'] = {
            'workers_used': self.max_workers,
            'processing_method': 'parallel'
        }
        
        # Add transcript resolution summary
        if processed_input.missense_variants:
            variants_with_transcripts = sum(1 for v in processed_input.missense_variants if v.get('transcript_id'))
            variants_with_gencode_ids = sum(1 for v in processed_input.missense_variants if v.get('gencode_id'))
            output_data['transcript_resolution_summary'] = {
                'total_variants': len(processed_input.missense_variants),
                'variants_with_transcript_ids': variants_with_transcripts,
                'transcript_resolution_rate': variants_with_transcripts / len(processed_input.missense_variants),
                'variants_with_gencode_ids': variants_with_gencode_ids,
                'gencode_id_resolution_rate': variants_with_gencode_ids / len(processed_input.missense_variants),
                'processing_method': 'parallel'
            }
        
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2, default=str)
        
        self.logger.info(f"Processed input exported to {output_path}")

def main():
    """Example usage of the enhanced input processor"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Enhanced Input Processor for ITS4.2')
    parser.add_argument('--pdb', help='Path to PDB file')
    parser.add_argument('--variants', help='Path to variants JSON file')
    parser.add_argument('--clinical', help='Path to clinical data JSON file')
    parser.add_argument('--dna', help='DNA sequence')
    parser.add_argument('--rna', help='RNA sequence')
    parser.add_argument('--protein', help='Protein sequence')
    parser.add_argument('--output', default='processed_input.json', help='Output file path')
    parser.add_argument('--workers', type=int, default=8, help='Number of parallel workers')
    
    args = parser.parse_args()
    
    # Build inputs dictionary from arguments
    inputs = {}
    
    if args.pdb:
        inputs['pdb_file'] = args.pdb
    if args.variants:
        inputs['missense_variants'] = args.variants
    if args.clinical:
        inputs['clinical_data'] = args.clinical
    if args.dna:
        inputs['dna_sequence'] = args.dna
    if args.rna:
        inputs['rna_sequence'] = args.rna
    if args.protein:
        inputs['protein_sequence'] = args.protein
    
    # Check if any inputs were provided
    if not inputs:
        logger.warning("No inputs provided. Please provide at least one of the following:")
        logger.warning("  --pdb: Path to PDB file")
        logger.warning("  --variants: Path to variants JSON file")
        logger.warning("  --clinical: Path to clinical data JSON file")
        logger.warning("  --dna: DNA sequence")
        logger.warning("  --rna: RNA sequence")
        logger.warning("  --protein: Protein sequence")
        return
    
    # Initialize processor
    processor = EnhancedInputProcessor()
    processor.max_workers = 8
    if args.workers:
        processor.max_workers = args.workers
    
    # Process inputs
    processed_input = processor.process_inputs(inputs)
    
    # Get summary
    summary = processor.get_input_summary(processed_input)
    logger.info("Input Summary:")
    logger.info(json.dumps(summary, indent=2))
    
    # Export if valid
    if processed_input.validation_result and processed_input.validation_result.is_valid:
        processor.export_processed_input(processed_input, args.output)
        logger.info(f"Processed input exported successfully to {args.output}")
    else:
        logger.error("Input validation failed:")
        if processed_input.validation_result and processed_input.validation_result.errors:
            for error in processed_input.validation_result.errors:
                logger.error(f"  - {error}")
        else:
            logger.error("  - Unknown validation error or validation result missing.")

if __name__ == "__main__":
    main() 