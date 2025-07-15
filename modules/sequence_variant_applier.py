#!/usr/bin/env python3
"""
Sequence Variant Applier

A module for applying genetic variants to reference sequences using BioPython.
Supports DNA, RNA, and protein sequences with various variant types.

Example usage:
    applier = SequenceVariantApplier()
    results = applier.apply_variants_to_sequences("variant_impact_results.json")
"""

import os
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
from dataclasses import dataclass, field
from datetime import datetime
from collections import defaultdict
from Bio.Data.IUPACData import protein_letters_3to1

# Configure logging
from .logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)

# BioPython imports
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqUtils import nt_search
    from Bio.Data import CodonTable
    BIOPYTHON_AVAILABLE = True
except ImportError:
    logger.warning("BioPython not available. Install with: pip install biopython")
    BIOPYTHON_AVAILABLE = False


@dataclass
class VariantApplication:
    """Result of applying a variant to a sequence"""
    variant_id: str
    original_sequence: str
    mutated_sequence: str
    variant_type: str
    position: int
    reference_allele: str
    alternate_allele: str
    sequence_type: str  # 'dna', 'rna', 'protein'
    validation_status: str
    warnings: List[str] = field(default_factory=list)
    processing_timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


@dataclass
class SequenceVariantResults:
    """Complete results of applying variants to sequences"""
    input_file: str
    dna_sequences: Dict[str, VariantApplication] = field(default_factory=dict)
    rna_sequences: Dict[str, VariantApplication] = field(default_factory=dict)
    protein_sequences: Dict[str, VariantApplication] = field(default_factory=dict)
    combined_sequences: Dict[str, str] = field(default_factory=dict)
    processing_timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    summary: Dict[str, Any] = field(default_factory=dict)


class SequenceVariantApplier:
    """Applies genetic variants to reference sequences using BioPython"""
    
    def __init__(self):
        """Initialize the sequence variant applier"""
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython is required but not available")
        
        # Standard codon table for translation
        self.codon_table = CodonTable.standard_dna_table
        
        # Sequence type detection patterns
        self.dna_pattern = re.compile(r'^[ATCGN]+$', re.IGNORECASE)
        self.rna_pattern = re.compile(r'^[AUCGN]+$', re.IGNORECASE)
        self.protein_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY*]+$', re.IGNORECASE)
        
        logger.info("Sequence Variant Applier initialized")
    
    def apply_variants_to_sequences(self, input_file: str) -> SequenceVariantResults:
        """
        Apply all variants from the input file to reference sequences
        
        Args:
            input_file: Path to variant_impact_results.json file
            
        Returns:
            SequenceVariantResults with all applied variants
        """
        # Load input data
        with open(input_file, 'r') as f:
            input_data = json.load(f)
        
        variants = input_data.get('missense_variants', [])
        dna_sequence = input_data.get('dna_sequence')
        rna_sequence = input_data.get('rna_sequence')
        protein_sequence = input_data.get('protein_sequence')
        
        logger.info(f"Applying {len(variants)} variants to sequences")
        
        # Handle file paths for sequences
        if dna_sequence and os.path.exists(dna_sequence):
            logger.info(f"Loading DNA sequence from file: {dna_sequence}")
            try:
                with open(dna_sequence, 'r') as f:
                    dna_sequence = f.read().strip()
                logger.info(f"Loaded DNA sequence of length: {len(dna_sequence)}")
            except Exception as e:
                logger.error(f"Error loading DNA sequence from {dna_sequence}: {e}")
                dna_sequence = None
        
        if rna_sequence and os.path.exists(rna_sequence):
            logger.info(f"Loading RNA sequence from file: {rna_sequence}")
            try:
                with open(rna_sequence, 'r') as f:
                    rna_sequence = f.read().strip()
                logger.info(f"Loaded RNA sequence of length: {len(rna_sequence)}")
            except Exception as e:
                logger.error(f"Error loading RNA sequence from {rna_sequence}: {e}")
                rna_sequence = None
        
        if protein_sequence and os.path.exists(protein_sequence):
            logger.info(f"Loading protein sequence from file: {protein_sequence}")
            try:
                with open(protein_sequence, 'r') as f:
                    protein_sequence = f.read().strip()
                logger.info(f"Loaded protein sequence of length: {len(protein_sequence)}")
            except Exception as e:
                logger.error(f"Error loading protein sequence from {protein_sequence}: {e}")
                protein_sequence = None
        
        # Initialize results
        results = SequenceVariantResults(input_file=input_file)
        
        # Generate combined sequences (all variants applied)
        results.combined_sequences = self._generate_combined_sequences(
            variants, dna_sequence, rna_sequence, protein_sequence
        )
        
        # Generate summary
        results.summary = self._generate_summary(results, variants)
        
        return results
    
    def _apply_variants_to_sequence(self, variants: List[Dict], 
                                   reference_sequence: str, 
                                   sequence_type: str) -> Dict[str, VariantApplication]:
        """
        Apply variants to a single reference sequence
        
        Args:
            variants: List of variant dictionaries
            reference_sequence: Reference sequence to apply variants to
            sequence_type: Type of sequence ('dna', 'rna', 'protein')
            
        Returns:
            Dictionary mapping variant_id to VariantApplication result
        """
        results = {}
        
        for variant in variants:
            variant_id = variant['id']
            logger.info(f"Applying variant {variant_id} to {sequence_type} sequence")
            
            try:
                # Apply single variant
                result = self._apply_single_variant(
                    variant, reference_sequence, sequence_type
                )
                results[variant_id] = result
                
            except Exception as e:
                logger.error(f"Error applying variant {variant_id}: {e}")
                # Create error result
                results[variant_id] = VariantApplication(
                    variant_id=variant_id,
                    original_sequence=reference_sequence,
                    mutated_sequence=reference_sequence,  # Keep original on error
                    variant_type='error',
                    position=variant.get('position', 0),
                    reference_allele=variant.get('reference', ''),
                    alternate_allele=variant.get('alternate', ''),
                    sequence_type=sequence_type,
                    validation_status='error',
                    warnings=[f"Failed to apply variant: {str(e)}"]
                )
        
        return results
    
    def _apply_single_variant(self, variant: Dict, reference_sequence: str, 
                             sequence_type: str) -> VariantApplication:
        """
        Apply a single variant to a reference sequence
        
        Args:
            variant: Variant dictionary
            reference_sequence: Reference sequence
            sequence_type: Type of sequence
            
        Returns:
            VariantApplication result
        """
        variant_id = variant['id']
        position = variant.get('position', 0)
        reference_allele = variant.get('reference', '')
        alternate_allele = variant.get('alternate', '')
        
        # Handle genomic coordinates vs local coordinates
        local_position = self._convert_genomic_to_local_position(
            position, reference_sequence, variant_id
        )
        
        # Log what we're looking for vs what we found
        if local_position >= 0 and local_position < len(reference_sequence):
            found_sequence = reference_sequence[local_position:local_position + len(reference_allele)] if reference_allele else reference_sequence[local_position]
            logger.info(f"Variant {variant_id}: expected '{reference_allele}' at position {position}, found '{found_sequence}'")
        
        # Validate variant
        validation_result = self._validate_variant(
            variant, reference_sequence, sequence_type, local_position
        )
        
        if not validation_result['is_valid']:
            return VariantApplication(
                variant_id=variant_id,
                original_sequence=reference_sequence,
                mutated_sequence=reference_sequence,
                variant_type='invalid',
                position=position,
                reference_allele=reference_allele,
                alternate_allele=alternate_allele,
                sequence_type=sequence_type,
                validation_status='invalid',
                warnings=validation_result['warnings']
            )
        
        # Apply the variant using local position
        mutated_sequence = self._apply_variant_to_sequence(
            reference_sequence, local_position, reference_allele, 
            alternate_allele, sequence_type
        )
        
        # Determine variant type
        variant_type = self._determine_variant_type(reference_allele, alternate_allele)
        
        return VariantApplication(
            variant_id=variant_id,
            original_sequence=reference_sequence,
            mutated_sequence=mutated_sequence,
            variant_type=variant_type,
            position=position,
            reference_allele=reference_allele,
            alternate_allele=alternate_allele,
            sequence_type=sequence_type,
            validation_status='valid',
            warnings=validation_result['warnings']
        )
    
    def _convert_genomic_to_local_position(self, genomic_position: int, 
                                         reference_sequence: str, 
                                         variant_id: str) -> int:
        """
        Convert genomic position to local sequence position (1-based to 0-based)
        
        Args:
            genomic_position: Genomic coordinate (1-based)
            reference_sequence: Reference sequence
            variant_id: Variant ID for logging
            
        Returns:
            Local position (0-based)
        """
        sequence_length = len(reference_sequence)
        
        # Convert 1-based to 0-based
        local_position = genomic_position - 1
        
        if local_position >= sequence_length:
            if local_position > sequence_length * 10:
                local_position = local_position % sequence_length
                logger.warning(f"Genomic position {genomic_position} for variant {variant_id} "
                             f"is much larger than sequence length {sequence_length}. "
                             f"Using position {local_position} (modulo)")
            else:
                logger.warning(f"Genomic position {genomic_position} for variant {variant_id} "
                             f"exceeds sequence length {sequence_length}. "
                             f"Using actual position {local_position}")
        return local_position
    
    def _validate_variant(self, variant: Dict, reference_sequence: str, 
                         sequence_type: str, local_position: int) -> Dict[str, Any]:
        """
        Validate a variant for application to a sequence
        
        Args:
            variant: Variant dictionary
            reference_sequence: Reference sequence
            sequence_type: Type of sequence
            local_position: Local position in sequence
            
        Returns:
            Validation result dictionary
        """
        warnings = []
        is_valid = True
        
        reference_allele = variant.get('reference', '')
        alternate_allele = variant.get('alternate', '')
        
        # Check position bounds
        if local_position < 0 or local_position >= len(reference_sequence):
            warnings.append(f"Local position {local_position} out of bounds for sequence length {len(reference_sequence)}")
            is_valid = False
        
        # Check reference allele matches
        if reference_allele and reference_allele != '-':
            if local_position + len(reference_allele) > len(reference_sequence):
                warnings.append(f"Reference allele extends beyond sequence bounds")
                is_valid = False
            else:
                actual_reference = reference_sequence[local_position:local_position + len(reference_allele)]
                if actual_reference.upper() != reference_allele.upper():
                    warnings.append(f"Reference allele mismatch: expected '{reference_allele}', found '{actual_reference}'")
                    is_valid = False
        
        # Validate allele format for sequence type
        if sequence_type == 'dna':
            if not self._is_valid_dna_sequence(reference_allele) or not self._is_valid_dna_sequence(alternate_allele):
                warnings.append("Invalid DNA sequence in alleles")
                is_valid = False
        elif sequence_type == 'rna':
            if not self._is_valid_rna_sequence(reference_allele) or not self._is_valid_rna_sequence(alternate_allele):
                warnings.append("Invalid RNA sequence in alleles")
                is_valid = False
        elif sequence_type == 'protein':
            if not self._is_valid_protein_sequence(reference_allele) or not self._is_valid_protein_sequence(alternate_allele):
                warnings.append("Invalid protein sequence in alleles")
                is_valid = False
        
        return {
            'is_valid': is_valid,
            'warnings': warnings
        }
    
    def _apply_variant_to_sequence(self, reference_sequence: str, position: int,
                                  reference_allele: str, alternate_allele: str,
                                  sequence_type: str) -> str:
        """
        Apply a variant to a reference sequence
        
        Args:
            reference_sequence: Reference sequence
            position: Position of the variant (0-based)
            reference_allele: Reference allele
            alternate_allele: Alternate allele
            sequence_type: Type of sequence
            
        Returns:
            Mutated sequence
        """
        # Convert to BioPython Seq object
        seq = Seq(reference_sequence)
        
        if reference_allele == '-' or not reference_allele:
            # Insertion
            if alternate_allele != '-':
                seq = seq[:position] + Seq(alternate_allele) + seq[position:]
        elif alternate_allele == '-':
            # Deletion
            if position + len(reference_allele) <= len(seq):
                seq = seq[:position] + seq[position + len(reference_allele):]
        else:
            # Substitution
            if position + len(reference_allele) <= len(seq):
                seq = seq[:position] + Seq(alternate_allele) + seq[position + len(reference_allele):]
        
        return str(seq)
    
    def _determine_variant_type(self, reference_allele: str, alternate_allele: str) -> str:
        """
        Determine the type of variant
        
        Args:
            reference_allele: Reference allele
            alternate_allele: Alternate allele
            
        Returns:
            Variant type string
        """
        if reference_allele == '-' or not reference_allele:
            return 'insertion'
        elif alternate_allele == '-':
            return 'deletion'
        elif len(reference_allele) == len(alternate_allele):
            return 'substitution'
        else:
            return 'complex'
    
    def _generate_combined_sequences(self, variants: List[Dict], 
                                   dna_sequence: Optional[str],
                                   rna_sequence: Optional[str],
                                   protein_sequence: Optional[str]) -> Dict[str, str]:
        """
        Generate sequences with all variants applied
        
        Args:
            variants: List of variants
            dna_sequence: DNA reference sequence
            rna_sequence: RNA reference sequence
            protein_sequence: Protein reference sequence
            
        Returns:
            Dictionary of combined sequences
        """
        combined_sequences = {}
        
        # Apply all variants to DNA sequence
        if dna_sequence:
            combined_dna = dna_sequence
            for variant in variants:
                position = variant.get('position', 0)
                reference_allele = variant.get('reference', '')
                alternate_allele = variant.get('alternate', '')
                
                # Convert to local position
                local_position = self._convert_genomic_to_local_position(
                    position, combined_dna, variant['id']
                )
                
                if self._validate_variant(variant, combined_dna, 'dna', local_position)['is_valid']:
                    combined_dna = self._apply_variant_to_sequence(
                        combined_dna, local_position, reference_allele, alternate_allele, 'dna'
                    )
            combined_sequences['dna'] = combined_dna
        
        # Apply all variants to RNA sequence
        if rna_sequence:
            combined_rna = rna_sequence
            for variant in variants:
                position = variant.get('position', 0)
                reference_allele = variant.get('reference', '')
                alternate_allele = variant.get('alternate', '')
                
                # Convert to local position
                local_position = self._convert_genomic_to_local_position(
                    position, combined_rna, variant['id']
                )
                
                if self._validate_variant(variant, combined_rna, 'rna', local_position)['is_valid']:
                    combined_rna = self._apply_variant_to_sequence(
                        combined_rna, local_position, reference_allele, alternate_allele, 'rna'
                    )
            combined_sequences['rna'] = combined_rna
        
        # Apply all variants to protein sequence
        if protein_sequence:
            combined_protein = protein_sequence
            for variant in variants:
                position = variant.get('position', 0)
                reference_allele = variant.get('reference', '')
                alternate_allele = variant.get('alternate', '')
                
                # Convert to local position
                local_position = self._convert_genomic_to_local_position(
                    position, combined_protein, variant['id']
                )
                
                if self._validate_variant(variant, combined_protein, 'protein', local_position)['is_valid']:
                    combined_protein = self._apply_variant_to_sequence(
                        combined_protein, local_position, reference_allele, alternate_allele, 'protein'
                    )
            combined_sequences['protein'] = combined_protein
        
        return combined_sequences
    
    def _generate_summary(self, results: SequenceVariantResults, 
                         variants: List[Dict]) -> Dict[str, Any]:
        """
        Generate summary statistics for the variant application
        
        Args:
            results: SequenceVariantResults object
            variants: List of variants
            
        Returns:
            Summary dictionary
        """
        total_variants = len(variants)
        variant_types = defaultdict(int)
        
        # Count variant types
        for variant in variants:
            reference_allele = variant.get('reference', '')
            alternate_allele = variant.get('alternate', '')
            variant_type = self._determine_variant_type(reference_allele, alternate_allele)
            variant_types[variant_type] += 1
        
        # Count sequences processed
        sequences_processed = {
            'dna': 1 if 'dna' in results.combined_sequences else 0,
            'rna': 1 if 'rna' in results.combined_sequences else 0,
            'protein': 1 if 'protein' in results.combined_sequences else 0
        }
        
        return {
            'total_variants': total_variants,
            'variant_type_distribution': dict(variant_types),
            'sequences_processed': sequences_processed,
            'combined_sequences_generated': len(results.combined_sequences)
        }
    
    def _is_valid_dna_sequence(self, sequence: str) -> bool:
        """Check if sequence is valid DNA"""
        if not sequence or sequence == '-':
            return True
        return bool(self.dna_pattern.match(sequence))
    
    def _is_valid_rna_sequence(self, sequence: str) -> bool:
        """Check if sequence is valid RNA"""
        if not sequence or sequence == '-':
            return True
        return bool(self.rna_pattern.match(sequence))
    
    def _is_valid_protein_sequence(self, sequence: str) -> bool:
        """Check if sequence is valid protein"""
        if not sequence or sequence == '-':
            return True
        return bool(self.protein_pattern.match(sequence))
    
    def save_results(self, results: SequenceVariantResults, output_file: str):
        """
        Save results to JSON file (only combined sequence and essential metadata)
        
        Args:
            results: SequenceVariantResults object
            output_file: Output file path
        """
        # Convert dataclass to dictionary with only essential data
        results_dict = {
            'input_file': results.input_file,
            'processing_timestamp': results.processing_timestamp,
            'summary': results.summary,
            # Only save combined sequences (all variants applied together)
            'combined_sequences': results.combined_sequences
        }
        
        with open(output_file, 'w') as f:
            json.dump(results_dict, f, indent=2)
        
        logger.info(f"Results saved to {output_file}")

    def read_clean_protein_sequence(self, file_path: str) -> str:
        """
        Read a protein sequence from a text file, remove non-alphabetic characters, and join lines.
        Args:
            file_path: Path to the protein sequence text file
        Returns:
            Cleaned protein sequence string
        """
        with open(file_path, 'r') as f:
            raw = f.read()
        # Remove all non-alphabetic characters (keep only A-Z, a-z)
        cleaned = ''.join([c for c in raw if c.isalpha()])
        return cleaned


def main():
    """Main CLI interface for sequence variant applier"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Apply variants to reference sequences")
    parser.add_argument("input_file", nargs="?", help="Input JSON file with variants", default="processed_input.json")
    parser.add_argument("-o", "--output", help="Output JSON file", 
                       default="sequence_variant_results.json")
    parser.add_argument("--protein-seq-file", help="Path to protein sequence text file", 
                       default="examples/ex_protein_seq.txt")
    
    args = parser.parse_args()
    
    # Initialize applier
    applier = SequenceVariantApplier()
    
    # Read and clean protein sequence
    protein_sequence = applier.read_clean_protein_sequence(args.protein_seq_file)
    logger.info(f"Loaded protein sequence of length: {len(protein_sequence)}")
    
    # Load input data (variants)
    with open(args.input_file, 'r') as f:
        input_data = json.load(f)
    variants = input_data.get('missense_variants', [])
    logger.info(f"Loaded {len(variants)} variants")
    
    # Create a modified input data with the protein sequence
    modified_input_data = input_data.copy()
    modified_input_data['protein_sequence'] = protein_sequence
    
    # Save modified input data to a temporary file
    temp_input_file = "temp_input.json"
    with open(temp_input_file, 'w') as f:
        json.dump(modified_input_data, f)
    
    try:
        # Use the existing apply_variants_to_sequences method
        results = applier.apply_variants_to_sequences(temp_input_file)
        
        # Save results
        applier.save_results(results, args.output)
        
        # Print summary
        logger.info(f"Processed {results.summary['total_variants']} variants")
        logger.info(f"Generated {results.summary['combined_sequences_generated']} combined sequences")
        logger.info(f"Results saved to {args.output}")
        
    finally:
        # Clean up temporary file
        if os.path.exists(temp_input_file):
            os.remove(temp_input_file)


if __name__ == "__main__":
    main() 