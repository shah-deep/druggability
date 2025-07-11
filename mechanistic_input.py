#!/usr/bin/env python3
"""
Enhanced Input Router for ITS4.2 Mechanistic Coherence Module
Extends the existing pocket detection pipeline to handle multi-modal inputs
"""

import os
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from dataclasses import dataclass
import logging

@dataclass
class MechanisticInput:
    """Data structure for mechanistic coherence inputs"""
    # Existing pocket data
    pocket_data: Optional[Dict] = None
    pdb_file: Optional[str] = None
    
    # New input types for ITS4.2
    dna_sequence: Optional[str] = None
    rna_sequence: Optional[str] = None
    protein_sequence: Optional[str] = None
    structural_constraints: Optional[Dict] = None
    single_cell_data: Optional[str] = None  # Path to h5ad file
    missense_variants: Optional[List[Dict]] = None
    clinical_data: Optional[Dict] = None

class ITS42InputRouter:
    """Enhanced input router for mechanistic coherence analysis"""
    
    def __init__(self, pocket_detector=None, feature_extractor=None):
        self.pocket_detector = pocket_detector
        self.feature_extractor = feature_extractor
        self.logger = logging.getLogger(__name__)
        
    def route_inputs(self, inputs: MechanisticInput) -> Dict[str, List[str]]:
        """Route inputs to appropriate model categories"""
        routing_map = {
            'genotype_models': [],
            'simulation_models': [],
            'structure_models': [],
            'biomarker_models': [],
            'metabolic_models': []
        }
        
        # Route based on available data types
        if inputs.dna_sequence or inputs.rna_sequence:
            routing_map['genotype_models'].extend(['Enformer', 'BigRNA', 'GeneGenie'])
            
        if inputs.missense_variants:
            routing_map['genotype_models'].append('AlphaMissense')
            
        if inputs.single_cell_data:
            routing_map['simulation_models'].extend(['scGPT', 'scFoundation'])
            
        if inputs.clinical_data:
            routing_map['simulation_models'].append('OCTO')
            routing_map['biomarker_models'].append('ArteraAI')
            
        if inputs.protein_sequence or inputs.pdb_file:
            routing_map['structure_models'].extend(['ProteinMPNN', 'RFDiffusion', 'ESM-3'])
            
        if inputs.rna_sequence:
            routing_map['structure_models'].append('ATOM-1')
            
        if inputs.dna_sequence:
            routing_map['metabolic_models'].append('COBRApy')
            
        return routing_map
    
    def prepare_model_inputs(self, inputs: MechanisticInput) -> Dict[str, Dict]:
        """Prepare inputs for each model category"""
        prepared_inputs = {}
        
        # Genotype model inputs
        if inputs.dna_sequence:
            prepared_inputs['enformer'] = {'sequence': inputs.dna_sequence}
            prepared_inputs['genegenie'] = {'sequence': inputs.dna_sequence}
            
        if inputs.rna_sequence:
            prepared_inputs['bigrna'] = {'sequence': inputs.rna_sequence}
            
        if inputs.missense_variants:
            prepared_inputs['alphamissense'] = {'variants': inputs.missense_variants}
            
        # Structure model inputs
        if inputs.protein_sequence:
            prepared_inputs['proteinmpnn'] = {'sequence': inputs.protein_sequence}
            prepared_inputs['esm3'] = {'sequence': inputs.protein_sequence}
            
        if inputs.pdb_file:
            prepared_inputs['rfdiffusion'] = {'structure_file': inputs.pdb_file}
            
        return prepared_inputs
