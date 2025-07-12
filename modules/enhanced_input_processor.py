#!/usr/bin/env python3
"""
Enhanced Input Processor for ITS4.2 Mechanistic Coherence Module
Handles multi-modal input validation, preprocessing, and routing for mechanistic coherence analysis
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any
from dataclasses import dataclass, field
import logging
import hashlib
from datetime import datetime
import warnings

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
        
    def _get_default_config(self) -> Dict:
        """Get default configuration"""
        return {
            'max_sequence_length': 1000000,
            'allowed_file_types': ['.pdb', '.json', '.h5ad', '.csv', '.txt'],
            'required_fields': ['pdb_file'],  # Core structural requirement
            'optional_fields': [
                'dna_sequence', 'rna_sequence', 'protein_sequence',
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
                'required': True,
                'file_exists': True,
                'file_extension': '.pdb',
                'max_size_mb': 100
            },
            'dna_sequence': {
                'required': False,
                'pattern': r'^[ATGCatgc\s]+$',
                'min_length': 10,
                'max_length': 1000000
            },
            'rna_sequence': {
                'required': False,
                'pattern': r'^[AUGCaugc\s]+$',
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
            processed.dna_sequence = self._normalize_sequence(inputs['dna_sequence'])
            
        if 'rna_sequence' in inputs:
            processed.rna_sequence = self._normalize_sequence(inputs['rna_sequence'])
            
        if 'protein_sequence' in inputs:
            processed.protein_sequence = self._normalize_sequence(inputs['protein_sequence'])
            
        return processed
    
    def _process_variant_inputs(self, inputs: Dict[str, Any], processed: ProcessedInput) -> ProcessedInput:
        """Process variant inputs"""
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
            
            processed.missense_variants = self._normalize_variants(variants)
            
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
                'transcript': str(variant.get('transcript', ''))
            }
            
            normalized.append(normalized_variant)
        
        return normalized
    
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
            'validation_passed': processed_input.validation_result.is_valid if processed_input.validation_result else False
        }
        
        # Add counts
        if processed_input.missense_variants:
            summary['variant_count'] = len(processed_input.missense_variants)
        
        if processed_input.dna_sequence:
            summary['dna_length'] = len(processed_input.dna_sequence)
        
        if processed_input.rna_sequence:
            summary['rna_length'] = len(processed_input.rna_sequence)
        
        if processed_input.protein_sequence:
            summary['protein_length'] = len(processed_input.protein_sequence)
        
        return summary
    
    def export_processed_input(self, processed_input: ProcessedInput, output_path: str) -> None:
        """Export processed input to file"""
        output_data = {
            'pdb_file': processed_input.pdb_file,
            'pocket_data': processed_input.pocket_data,
            'dna_sequence': processed_input.dna_sequence,
            'rna_sequence': processed_input.rna_sequence,
            'protein_sequence': processed_input.protein_sequence,
            'missense_variants': processed_input.missense_variants,
            'clinical_data': processed_input.clinical_data,
            'structural_constraints': processed_input.structural_constraints,
            'single_cell_data': processed_input.single_cell_data,
            'metabolic_data': processed_input.metabolic_data,
            'input_hash': processed_input.input_hash,
            'processing_timestamp': processed_input.processing_timestamp,
            'validation_result': {
                'is_valid': processed_input.validation_result.is_valid,
                'errors': processed_input.validation_result.errors,
                'warnings': processed_input.validation_result.warnings
            } if processed_input.validation_result else None
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
    
    # Initialize processor
    processor = EnhancedInputProcessor()
    
    # Process inputs
    processed_input = processor.process_inputs(inputs)
    
    # Get summary
    summary = processor.get_input_summary(processed_input)
    print("Input Summary:")
    print(json.dumps(summary, indent=2))
    
    # Export if valid
    if processed_input.validation_result and processed_input.validation_result.is_valid:
        processor.export_processed_input(processed_input, args.output)
        print(f"Processed input exported successfully to {args.output}")
    else:
        print("Input validation failed:")
        if processed_input.validation_result and processed_input.validation_result.errors:
            for error in processed_input.validation_result.errors:
                print(f"  - {error}")
        else:
            print("  - Unknown validation error or validation result missing.")

if __name__ == "__main__":
    main() 