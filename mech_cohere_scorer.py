"""
Mechanistic Coherence Scoring System
Integrates all model outputs into coherence scores
"""

import numpy as np
from typing import Dict, List
import logging

class MechanisticCoherenceScorer:
    """Constructs mechanistic coherence scores from model outputs"""
    
    def __init__(self, weights=None):
        # Default weights for each category
        self.weights = weights or {
            'genotype': 0.25,
            'simulation': 0.25,
            'structure': 0.20,
            'biomarker': 0.15,
            'metabolic': 0.15
        }
        
        # Ensure weights sum to 1
        total_weight = sum(self.weights.values())
        self.weights = {k: v/total_weight for k, v in self.weights.items()}
    
    def calculate_category_scores(self, model_outputs: Dict) -> Dict:
        """Calculate scores for each mechanistic category"""
        category_scores = {}
        
        # Genotype score
        genotype_outputs = model_outputs.get('genotype', [])
        if genotype_outputs:
            genotype_score = np.mean([
                output.get('expression_score', 0.5) if 'expression_score' in output
                else output.get('splicing_score', 0.5) if 'splicing_score' in output  
                else 1 - output.get('overall_pathogenicity', 0.5) if 'overall_pathogenicity' in output
                else output.get('regulatory_score', 0.5)
                for output in genotype_outputs
            ])
        else:
            genotype_score = 0.5
        
        category_scores['genotype'] = genotype_score
        
        # Simulation score
        simulation_outputs = model_outputs.get('simulation', [])
        if simulation_outputs:
            simulation_score = np.mean([
                output.get('outcome_probability', 0.5) if 'outcome_probability' in output
                else output.get('cell_viability', 0.5) if 'cell_viability' in output
                else output.get('drug_efficacy', 0.5) if 'drug_efficacy' in output
                else output.get('network_score', 0.5) if 'network_score' in output
                else output.get('pathway_score', 0.5)
                for output in simulation_outputs
            ])
        else:
            simulation_score = 0.5
            
        category_scores['simulation'] = simulation_score
        
        # Structure score
        structure_outputs = model_outputs.get('structure', [])
        if structure_outputs:
            structure_score = np.mean([
                output.get('confidence', 0.5) if 'confidence' in output
                else output.get('design_score', 0.5) if 'design_score' in output
                else output.get('structure_confidence', 0.5) if 'structure_confidence' in output
                else output.get('evolutionary_score', 0.5)
                for output in structure_outputs
            ])
        else:
            structure_score = 0.5
            
        category_scores['structure'] = structure_score
        
        # Biomarker score
        biomarker_outputs = model_outputs.get('biomarker', [])
        if biomarker_outputs:
            biomarker_score = np.mean([
                output.get('biomarker_score', 0.5) if 'biomarker_score' in output
                else output.get('overall_feasibility', 0.5)
                for output in biomarker_outputs
            ])
        else:
            biomarker_score = 0.5
            
        category_scores['biomarker'] = biomarker_score
        
        # Metabolic score
        metabolic_outputs = model_outputs.get('metabolic', [])
        if metabolic_outputs:
            metabolic_score = np.mean([
                output.get('pathway_resilience', 0.5) if 'pathway_resilience' in output
                else output.get('metabolic_efficiency', 0.5) if 'metabolic_efficiency' in output
                else output.get('metabolic_safety_score', 0.5)
                for output in metabolic_outputs
            ])
        else:
            metabolic_score = 0.5
            
        category_scores['metabolic'] = metabolic_score
        
        return category_scores
    
    def calculate_final_score(self, category_scores: Dict) -> float:
        """Calculate weighted final coherence score"""
        final_score = sum(
            category_scores.get(category, 0.5) * weight 
            for category, weight in self.weights.items()
        )
        return final_score
    
    def detect_conflicts(self, model_outputs: Dict) -> List[str]:
        """Detect conflicts between model predictions"""
        conflicts = []
        
        # Check for genotype conflicts
        genotype_outputs = model_outputs.get('genotype', [])
        enformer_output = next((o for o in genotype_outputs if o.get('model') == 'Enformer'), None)
        alphamissense_output = next((o for o in genotype_outputs if o.get('model') == 'AlphaMissense'), None)
        
        if enformer_output and alphamissense_output:
            expr_score = enformer_output.get('expression_score', 0.5)
            pathogenicity = alphamissense_output.get('overall_pathogenicity', 0.5)
            
            if expr_score > 0.6 and pathogenicity > 0.7:
                conflicts.append("Enformer vs AlphaMissense: Normal expression despite high pathogenicity")
        
        # Check for structure-function conflicts
        structure_outputs = model_outputs.get('structure', [])
        proteinmpnn_output = next((o for o in structure_outputs if o.get('model') == 'ProteinMPNN'), None)
        esm3_output = next((o for o in structure_outputs if o.get('model') == 'ESM-3'), None)
        
        if proteinmpnn_output and esm3_output:
            binding_affinity = proteinmpnn_output.get('binding_affinity_estimate', 0.5)
            evolutionary_score = esm3_output.get('evolutionary_score', 0.5)
            
            if binding_affinity > 0.8 and evolutionary_score < 0.3:
                conflicts.append("ProteinMPNN vs ESM-3: High binding affinity but low evolutionary plausibility")
        
        return conflicts
    
    def calculate_confidence(self, model_outputs: Dict, conflicts: List[str]) -> float:
        """Calculate overall confidence based on model agreement"""
        # Start with base confidence
        base_confidence = 0.8
        
        # Reduce confidence for each conflict
        conflict_penalty = len(conflicts) * 0.1
        
        # Calculate confidence based on model consistency
        all_confidences = []
        for category_outputs in model_outputs.values():
            for output in category_outputs:
                if 'confidence' in output:
                    all_confidences.append(output['confidence'])
        
        if all_confidences:
            avg_model_confidence = np.mean(all_confidences)
            model_agreement = 1.0 - np.std(all_confidences)  # Higher std = lower agreement
        else:
            avg_model_confidence = 0.7
            model_agreement = 0.8
        
        final_confidence = base_confidence * avg_model_confidence * model_agreement - conflict_penalty
        return max(0.1, min(1.0, final_confidence))