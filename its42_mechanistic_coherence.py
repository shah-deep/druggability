#!/usr/bin/env python3
"""
ITS4.2 Mechanistic Coherence Module
Building upon existing druggability pipeline from Task 1, also made by Joshua Robert
Implemented with ~10 of suggested models from instruction documentation

Author: Joshua Robert
Date: July 2025
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional

# Import existing modules
from pocket_detection import PocketDetector
from feature_extraction import FeatureExtractor

# Import new ITS4.2 modules
from mechanistic_input import ITS42InputRouter, MechanisticInput
from genotype_predictor import GenotypePredictor
from mechanistics_simulator import MechanisticSimulator
from struc_func_predictor import StructureFunctionPredictor
from biomarker_anchor_system import BiomarkerAnchoringSystem
from metab_cohere_analyzer import MetabolicCoherenceAnalyzer
from mech_cohere_scorer import MechanisticCoherenceScorer
from mech_cohere_output import MechanisticCoherenceOutput

class ITS42MechanisticCoherence:
    """
    Complete ITS4.2 Mechanistic Coherence Module
    Integrates all model categories and produces coherence scores
    """
    
    def __init__(self):
        # Initialize existing Task 1 components
        self.pocket_detector = PocketDetector()
        self.feature_extractor = FeatureExtractor()
        
        # Initialize ITS4.2 components
        self.input_router = ITS42InputRouter(self.pocket_detector, self.feature_extractor)
        self.genotype_predictor = GenotypePredictor()
        self.mechanistic_simulator = MechanisticSimulator()
        self.structure_predictor = StructureFunctionPredictor()
        self.biomarker_system = BiomarkerAnchoringSystem()
        self.metabolic_analyzer = MetabolicCoherenceAnalyzer()
        self.scorer = MechanisticCoherenceScorer()
        self.output_system = MechanisticCoherenceOutput(self.scorer)
        
    def analyze_mechanistic_coherence(self, inputs: MechanisticInput, molecule_id: str = None) -> Dict:
        """
        Complete mechanistic coherence analysis pipeline
        """
        print(f"Starting mechanistic coherence analysis for {molecule_id or 'unknown molecule'}")
        
        # Phase 1: Route inputs
        routing_map = self.input_router.route_inputs(inputs)
        prepared_inputs = self.input_router.prepare_model_inputs(inputs)
        
        print(f"Routing inputs to {sum(len(models) for models in routing_map.values())} models")
        
        # Phase 2: Run models by category
        model_outputs = {}
        
        # (A) Genotype → Phenotype Predictors
        if routing_map['genotype_models']:
            print("Running genotype prediction models...")
            genotype_outputs = []
            
            if inputs.dna_sequence and 'Enformer' in routing_map['genotype_models']:
                result = self.genotype_predictor.predict_expression(inputs.dna_sequence)
                genotype_outputs.append(result)
            
            if inputs.rna_sequence and 'BigRNA' in routing_map['genotype_models']:
                result = self.genotype_predictor.predict_splicing(inputs.rna_sequence)
                genotype_outputs.append(result)
            
            if inputs.missense_variants and 'AlphaMissense' in routing_map['genotype_models']:
                result = self.genotype_predictor.predict_variant_pathogenicity(inputs.missense_variants)
                genotype_outputs.append(result)
            
            if inputs.dna_sequence and 'GeneGenie' in routing_map['genotype_models']:
                result = self.genotype_predictor.predict_regulatory_activity(inputs.dna_sequence)
                genotype_outputs.append(result)
            
            model_outputs['genotype'] = genotype_outputs
        
        # (B) Mechanistic Simulations
        if routing_map['simulation_models']:
            print("Running mechanistic simulation models...")
            simulation_outputs = []
            
            if inputs.clinical_data and 'OCTO' in routing_map['simulation_models']:
                result = self.mechanistic_simulator.predict_patient_outcomes(inputs.clinical_data)
                simulation_outputs.append(result)
            
            if inputs.single_cell_data and 'scGPT' in routing_map['simulation_models']:
                result = self.mechanistic_simulator.predict_single_cell_response(inputs.single_cell_data)
                simulation_outputs.append(result)
            
            if inputs.single_cell_data and 'scFoundation' in routing_map['simulation_models']:
                # Mock expression data for demonstration
                expression_data = {'genes': ['TP53', 'BRCA1'], 'expression': [0.8, 0.6]}
                result = self.mechanistic_simulator.predict_transcriptomic_response(expression_data)
                simulation_outputs.append(result)
            
            model_outputs['simulation'] = simulation_outputs
        
        # (C) Protein/RNA Structure and Function
        if routing_map['structure_models']:
            print("Running structure and function models...")
            structure_outputs = []
            
            # Pass pocket data to structure predictor
            if inputs.pocket_data:
                self.structure_predictor.pocket_data = inputs.pocket_data
            
            if inputs.pdb_file and 'ProteinMPNN' in routing_map['structure_models']:
                result = self.structure_predictor.design_protein_sequence(inputs.pdb_file)
                structure_outputs.append(result)
            
            if inputs.pdb_file and 'RFDiffusion' in routing_map['structure_models']:
                constraints = {'binding_site': inputs.pocket_data.get('center', [0,0,0]) if inputs.pocket_data else [0,0,0]}
                result = self.structure_predictor.generate_protein_structure(constraints)
                structure_outputs.append(result)
            
            if inputs.protein_sequence and 'ESM-3' in routing_map['structure_models']:
                result = self.structure_predictor.assess_evolutionary_plausibility(inputs.protein_sequence)
                structure_outputs.append(result)
            
            if inputs.rna_sequence and 'ATOM-1' in routing_map['structure_models']:
                result = self.structure_predictor.predict_rna_structure(inputs.rna_sequence)
                structure_outputs.append(result)
            
            model_outputs['structure'] = structure_outputs
        
        # (D) Biomarker Anchoring
        if routing_map['biomarker_models']:
            print("Running biomarker anchoring models...")
            biomarker_outputs = []
            
            if inputs.clinical_data and 'ArteraAI' in routing_map['biomarker_models']:
                result = self.biomarker_system.predict_biomarkers(inputs.clinical_data)
                biomarker_outputs.append(result)
            
            # Mock target sites for OpenCRISPR-1
            if 'OpenCRISPR-1' in routing_map['biomarker_models']:
                target_sites = ['GGTCACGTACGTGCCAGG', 'AGTCGATCGATCGTAGCC']
                result = self.biomarker_system.assess_editing_feasibility(target_sites)
                biomarker_outputs.append(result)
            
            model_outputs['biomarker'] = biomarker_outputs
        
        # (E) Metabolic Pathway Coherence
        if routing_map['metabolic_models']:
            print("Running metabolic coherence analysis...")
            metabolic_outputs = []
            
            # Mock gene targets and drug targets
            gene_targets = ['TP53', 'BRCA1', 'EGFR']
            drug_targets = ['PI3K', 'mTOR', 'PARP1']
            
            knockout_result = self.metabolic_analyzer.simulate_gene_knockout(gene_targets)
            metabolic_outputs.append(knockout_result)
            
            drug_effect_result = self.metabolic_analyzer.predict_metabolic_drug_effects(drug_targets)
            metabolic_outputs.append(drug_effect_result)
            
            model_outputs['metabolic'] = metabolic_outputs
        
        # Phase 3: Score Construction
        print("Calculating mechanistic coherence scores...")
        
        # Phase 4: Output Formatting
        coherence_report = self.output_system.generate_coherence_report(
            model_outputs, molecule_id or "unknown"
        )
        
        print(f"Analysis complete. Final coherence score: {coherence_report['mechanistic_coherence_score']:.3f}")
        
        return coherence_report
    
    def run_complete_pipeline(self, pdb_file: str, **kwargs) -> Dict:
        """
        Run complete pipeline including Task 1 druggability + ITS4.2 mechanistic coherence
        """
        molecule_id = Path(pdb_file).stem
        print(f"Running complete pipeline for {molecule_id}")
        
        # Phase 1: Run existing Task 1 pipeline
        print("Running Task 1 druggability analysis...")
        
        # Detect pockets
        pockets, pocket_output_dir = self.pocket_detector.detect_pockets(pdb_file)
        
        # Extract features
        pocket_json_file = f"pockets_{molecule_id}.json"
        self.pocket_detector.save_pockets_json(pockets, pocket_json_file)
        
        with open(f"pocket_features_{molecule_id}.json") as f:
            features = json.load(f)
        
        # Phase 2: Run ITS4.2 mechanistic coherence
        print("Running ITS4.2 mechanistic coherence analysis...")
        
        # Prepare inputs
        mechanistic_inputs = MechanisticInput(
            pocket_data=pockets[0] if pockets else None,  # Use best pocket
            pdb_file=pdb_file,
            dna_sequence=kwargs.get('dna_sequence'),
            rna_sequence=kwargs.get('rna_sequence'), 
            protein_sequence=kwargs.get('protein_sequence'),
            structural_constraints=kwargs.get('structural_constraints'),
            single_cell_data=kwargs.get('single_cell_data'),
            missense_variants=kwargs.get('missense_variants'),
            clinical_data=kwargs.get('clinical_data')
        )
        
        # Run mechanistic coherence analysis
        coherence_report = self.analyze_mechanistic_coherence(mechanistic_inputs, molecule_id)
        
        # Phase 3: Integrate results
        print("Integrating druggability and mechanistic coherence results...")
        
        integrated_report = self.output_system.integrate_with_druggability(
            coherence_report, features
        )
        
        return integrated_report
    
    def save_results(self, results: Dict, output_file: str):
        """Save integrated results to JSON file"""
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Results saved to {output_file}")

def main():
    """Main CLI interface for ITS4.2"""
    parser = argparse.ArgumentParser(description="ITS4.2 Mechanistic Coherence Module")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("-o", "--output", help="Output JSON file", 
                       default="mechanistic_coherence_results.json")
    parser.add_argument("--dna-sequence", help="DNA sequence")
    parser.add_argument("--rna-sequence", help="RNA sequence") 
    parser.add_argument("--protein-sequence", help="Protein sequence")
    parser.add_argument("--single-cell-data", help="Single-cell data file (.h5ad)")
    parser.add_argument("--clinical-data", help="Clinical data JSON file")
    parser.add_argument("--missense-variants", help="Missense variants JSON file")
    
    args = parser.parse_args()
    
    # Initialize ITS4.2 system
    its42 = ITS42MechanisticCoherence()
    
    # Prepare optional inputs
    kwargs = {}
    if args.dna_sequence:
        kwargs['dna_sequence'] = args.dna_sequence
    if args.rna_sequence:
        kwargs['rna_sequence'] = args.rna_sequence
    if args.protein_sequence:
        kwargs['protein_sequence'] = args.protein_sequence
    if args.single_cell_data:
        kwargs['single_cell_data'] = args.single_cell_data
    if args.clinical_data:
        with open(args.clinical_data) as f:
            kwargs['clinical_data'] = json.load(f)
    if args.missense_variants:
        with open(args.missense_variants) as f:
            kwargs['missense_variants'] = json.load(f)
    
    try:
        # Run complete pipeline
        results = its42.run_complete_pipeline(args.pdb_file, **kwargs)
        
        # Save results
        its42.save_results(results, args.output)
        
        # Print summary
        coherence_score = results['compound_analysis']['mechanistic_coherence']['mechanistic_coherence_score']
        combined_score = results['combined_score']['combined_score']
        recommendation = results['recommendation']
        
        print(f"\n=== ITS4.2 Analysis Summary ===")
        print(f"Mechanistic Coherence Score: {coherence_score:.3f}")
        print(f"Combined Score: {combined_score:.3f}")  
        print(f"Recommendation: {recommendation}")
        print(f"Full results saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
