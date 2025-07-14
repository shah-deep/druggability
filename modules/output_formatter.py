#!/usr/bin/env python3
"""
Output Formatter for Pathogenicity Scoring Results

Converts comprehensive pathogenicity scoring results into standardized output format.
"""

import json
from typing import Dict, List, Any, Optional
from dataclasses import asdict
from comprehensive_pathogenicity_scorer import ComprehensiveScore

class OutputFormatter:
    """Format pathogenicity scoring results into standardized output"""
    
    def __init__(self):
        self.annotation_sources = {
            'AlphaMissense': 'alphamissense_score',
            'CADD': 'cadd_phred_score', 
            'VEP': 'vep_score',
            'ClinVar': 'clinvar_significance'
        }
    
    def format_variant_result(self, variant_id: str, score: ComprehensiveScore) -> Dict[str, Any]:
        """Format a single variant result into the specified output format"""
        
        # Calculate pathogenicity score (0-1 scale)
        pathogenicity_score = self._calculate_pathogenicity_score(score)
        
        # Calculate confidence based on available methods
        confidence = self._calculate_confidence(score)
        
        # Collect known annotations
        known_annotations = self._collect_annotations(score)
        
        return {
            "pathogenicity_score": pathogenicity_score,
            "confidence": confidence,
            "known_annotations": known_annotations
        }
    
    def format_all_results(self, results: Dict[str, ComprehensiveScore]) -> Dict[str, Dict[str, Any]]:
        """Format all variant results into the specified output format"""
        formatted_results = {}
        
        for variant_id, score in results.items():
            formatted_results[variant_id] = self.format_variant_result(variant_id, score)
        
        return formatted_results
    
    def _calculate_pathogenicity_score(self, score: ComprehensiveScore) -> float:
        """Calculate pathogenicity score (0-1 scale)"""
        if score.combined_score is not None:
            return score.combined_score
        
        # Fallback calculation if combined_score not available
        scores = []
        weights = []
        
        # AlphaMissense score
        if score.alphamissense_score is not None:
            scores.append(score.alphamissense_score)
            weights.append(0.4)
        
        # CADD score (normalize to 0-1)
        if score.cadd_phred_score is not None:
            # CADD scores are typically 0-50, normalize to 0-1
            normalized_cadd = min(score.cadd_phred_score / 50.0, 1.0)
            scores.append(normalized_cadd)
            weights.append(0.3)
        
        # VEP score
        if score.vep_score is not None:
            scores.append(score.vep_score)
            weights.append(0.2)
        
        # ClinVar significance (convert to numeric)
        if score.clinvar_significance:
            clinvar_score = self._clinvar_to_score(score.clinvar_significance)
            scores.append(clinvar_score)
            weights.append(0.1)
        
        if not scores:
            return 0.5  # Default neutral score
        
        # Weighted average
        total_weight = sum(weights)
        if total_weight > 0:
            weighted_score = sum(s * w for s, w in zip(scores, weights)) / total_weight
            return min(max(weighted_score, 0.0), 1.0)
        
        return 0.5
    
    def _calculate_confidence(self, score: ComprehensiveScore) -> float:
        """Calculate confidence score (0-1) based on available methods"""
        available_methods = len(score.methods_used)
        total_methods = len(self.annotation_sources)
        
        # Base confidence on method coverage
        method_coverage = available_methods / total_methods if total_methods > 0 else 0.0
        
        # Boost confidence if we have high-quality scores
        quality_boost = 0.0
        if score.alphamissense_score is not None:
            quality_boost += 0.2
        if score.cadd_phred_score is not None:
            quality_boost += 0.2
        if score.vep_score is not None:
            quality_boost += 0.1
        if score.clinvar_significance:
            quality_boost += 0.1
        
        confidence = min(method_coverage + quality_boost, 1.0)
        return round(confidence, 3)
    
    def _collect_annotations(self, score: ComprehensiveScore) -> List[Dict[str, Any]]:
        """Collect all known annotations for the variant"""
        annotations = []
        
        # AlphaMissense annotation
        if score.alphamissense_score is not None:
            annotations.append({
                "source": "AlphaMissense",
                "score": score.alphamissense_score,
                "type": "pathogenicity_prediction",
                "description": f"AlphaMissense pathogenicity score: {score.alphamissense_score:.3f}"
            })
        
        # CADD annotation
        if score.cadd_phred_score is not None:
            annotations.append({
                "source": "CADD",
                "score": score.cadd_phred_score,
                "type": "deleteriousness_prediction",
                "description": f"CADD Phred score: {score.cadd_phred_score:.2f}"
            })
        
        # VEP annotation
        if score.vep_score is not None:
            annotations.append({
                "source": "VEP",
                "score": score.vep_score,
                "type": "functional_impact",
                "description": f"VEP impact score: {score.vep_score:.3f}",
                "consequence": score.consequence,
                "impact": score.impact
            })
        
        # ClinVar annotation
        if score.clinvar_significance:
            annotations.append({
                "source": "ClinVar",
                "significance": score.clinvar_significance,
                "type": "clinical_significance",
                "description": f"ClinVar clinical significance: {score.clinvar_significance}"
            })
        
        # Variant type annotation
        if score.variant_type:
            annotations.append({
                "source": "Variant_Classification",
                "type": "variant_type",
                "description": f"Variant type: {score.variant_type}"
            })
        
        # Pathogenicity level annotation
        if score.pathogenicity_level:
            annotations.append({
                "source": "Combined_Analysis",
                "level": score.pathogenicity_level,
                "type": "pathogenicity_classification",
                "description": f"Pathogenicity level: {score.pathogenicity_level}"
            })
        
        return annotations
    
    def _clinvar_to_score(self, significance: str) -> float:
        """Convert ClinVar significance to numeric score"""
        significance_lower = significance.lower()
        
        if 'pathogenic' in significance_lower:
            return 0.9
        elif 'likely_pathogenic' in significance_lower:
            return 0.8
        elif 'uncertain' in significance_lower:
            return 0.5
        elif 'likely_benign' in significance_lower:
            return 0.2
        elif 'benign' in significance_lower:
            return 0.1
        else:
            return 0.5
    
    def export_formatted_results(self, results: Dict[str, ComprehensiveScore], output_path: str):
        """Export formatted results to JSON file"""
        formatted_results = self.format_all_results(results)
        
        # Add metadata
        output_data = {
            "metadata": {
                "format_version": "1.0",
                "description": "Pathogenicity scoring results in standardized format",
                "total_variants": len(formatted_results),
                "annotation_sources": list(self.annotation_sources.keys())
            },
            "results": formatted_results
        }
        
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        print(f"Formatted results exported to: {output_path}")
    
    def print_formatted_results(self, results: Dict[str, ComprehensiveScore]):
        """Print formatted results to console"""
        formatted_results = self.format_all_results(results)
        
        print("\n" + "="*60)
        print("FORMATTED PATHOGENICITY SCORING RESULTS")
        print("="*60)
        
        for variant_id, result in formatted_results.items():
            print(f"\nVariant: {variant_id}")
            print(f"  Pathogenicity Score: {result['pathogenicity_score']:.3f}")
            print(f"  Confidence: {result['confidence']:.3f}")
            print(f"  Annotations ({len(result['known_annotations'])}):")
            
            for annotation in result['known_annotations']:
                source = annotation.get('source', 'Unknown')
                description = annotation.get('description', 'No description')
                print(f"    - {source}: {description}")
        
        print("\n" + "="*60) 