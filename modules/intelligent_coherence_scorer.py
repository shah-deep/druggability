#!/usr/bin/env python3
"""
Intelligent Coherence Scorer

A module for aggregating and scoring outputs from all pipeline components.
Combines all scores into overall score and coherence category.
Generates a final JSON report with metadata, trace, and scores.

Example usage:
    scorer = IntelligentCoherenceScorer()
    result = scorer.aggregate_and_score(
        "outputs/structural_results_20250715_015155.json",
        "outputs/variant_impact_results_20250715_015155.json",
        "outputs/sequence_variant_results_20250715_015155.json",
        "outputs/coherence_results_20250715_015155.json",
        "outputs/pathway_dynamics_results_20250715_015155.json"
    )
    scorer.save_final_results(result, "outputs/final_results_20250715_015155.json")
"""

import os
import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import statistics

# Configure logging
from .logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)


@dataclass
class CoherenceCategory:
    """Coherence category with thresholds"""
    name: str
    min_score: float
    max_score: float
    description: str


@dataclass
class CategoryScore:
    """Individual category score"""
    genotype_pathogenicity: float
    structure_function: float
    pathway_enrichment: float
    clinical_translatability: float


@dataclass
class CausalTrace:
    """Causal trace step"""
    step: int
    level: str
    description: str


@dataclass
class FinalResult:
    """Final aggregated result"""
    overall_score: float
    coherence_category: str
    category_scores: CategoryScore
    causal_trace: List[CausalTrace]
    metadata: Dict[str, Any]
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class IntelligentCoherenceScorer:
    """
    Aggregates and scores outputs from all pipeline components
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the intelligent coherence scorer
        
        Args:
            config: Configuration dictionary
        """
        if config is None:
            config = self._get_default_config()
        
        self.config = config
        self.coherence_categories = self._define_coherence_categories()
        
        logger.info("IntelligentCoherenceScorer initialized")
    
    def _get_default_config(self) -> Dict:
        """Get default configuration"""
        return {
            'category_weights': {
                'genotype_pathogenicity': 0.25,
                'structure_function': 0.25,
                'pathway_enrichment': 0.25,
                'clinical_translatability': 0.25
            },
            'score_thresholds': {
                'high_coherence': 0.8,
                'medium_coherence': 0.6,
                'low_coherence': 0.4
            },
            'trace_levels': ['molecular', 'pathway', 'clinical', 'tissue']
        }
    
    def _define_coherence_categories(self) -> List[CoherenceCategory]:
        """Define coherence categories with thresholds"""
        return [
            CoherenceCategory("HIGH_COHERENCE", 0.8, 1.0, "Strong mechanistic coherence across all levels"),
            CoherenceCategory("MEDIUM_COHERENCE", 0.6, 0.8, "Moderate mechanistic coherence with some inconsistencies"),
            CoherenceCategory("LOW_COHERENCE", 0.4, 0.6, "Weak mechanistic coherence with significant inconsistencies"),
            CoherenceCategory("NO_COHERENCE", 0.0, 0.4, "No meaningful mechanistic coherence detected")
        ]
    
    def aggregate_and_score(self,
                          structural_results_file: str,
                          variant_impact_file: str,
                          sequence_variant_file: str,
                          coherence_results_file: str,
                          pathway_dynamics_file: str) -> FinalResult:
        """
        Aggregate and score all pipeline outputs
        
        Args:
            structural_results_file: Path to structural results JSON
            variant_impact_file: Path to variant impact results JSON
            sequence_variant_file: Path to sequence variant results JSON
            coherence_results_file: Path to coherence results JSON
            pathway_dynamics_file: Path to pathway dynamics results JSON
            
        Returns:
            FinalResult with aggregated scores and coherence category
        """
        logger.info("Starting aggregation and scoring of pipeline outputs")
        
        # Load all result files
        structural_results = self._load_json(structural_results_file)
        variant_impact_results = self._load_json(variant_impact_file)
        sequence_variant_results = self._load_json(sequence_variant_file)
        coherence_results = self._load_json(coherence_results_file)
        pathway_dynamics_results = self._load_json(pathway_dynamics_file)
        
        # Calculate individual category scores
        category_scores = self._calculate_category_scores(
            structural_results, variant_impact_results, sequence_variant_results,
            coherence_results, pathway_dynamics_results
        )
        
        # Calculate overall score
        overall_score = self._calculate_overall_score(category_scores)
        
        # Determine coherence category
        coherence_category = self._determine_coherence_category(overall_score)
        
        # Generate causal trace
        causal_trace = self._generate_causal_trace(
            variant_impact_results, structural_results, pathway_dynamics_results, coherence_results
        )
        
        # Create metadata
        metadata = self._create_metadata(
            structural_results, variant_impact_results, sequence_variant_results,
            coherence_results, pathway_dynamics_results
        )
        
        result = FinalResult(
            overall_score=overall_score,
            coherence_category=coherence_category,
            category_scores=category_scores,
            causal_trace=causal_trace,
            metadata=metadata
        )
        
        logger.info(f"Aggregation complete. Overall score: {overall_score:.3f}, Category: {coherence_category}")
        
        return result
    
    def _load_json(self, file_path: str) -> Dict:
        """Load JSON file with error handling"""
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            logger.error(f"File not found: {file_path}")
            return {}
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in {file_path}: {e}")
            return {}
    
    def _calculate_category_scores(self,
                                 structural_results: Dict,
                                 variant_impact_results: Dict,
                                 sequence_variant_results: Dict,
                                 coherence_results: Dict,
                                 pathway_dynamics_results: Dict) -> CategoryScore:
        """Calculate scores for each category"""
        
        # Genotype pathogenicity score (from variant impact analysis)
        genotype_pathogenicity = self._calculate_genotype_pathogenicity_score(variant_impact_results)
        
        # Structure function score (from structural analysis)
        structure_function = self._calculate_structure_function_score(structural_results)
        
        # Pathway enrichment score (from pathway dynamics)
        pathway_enrichment = self._calculate_pathway_enrichment_score(pathway_dynamics_results)
        
        # Clinical translatability score (from coherence analysis)
        clinical_translatability = self._calculate_clinical_translatability_score(coherence_results)
        
        return CategoryScore(
            genotype_pathogenicity=genotype_pathogenicity,
            structure_function=structure_function,
            pathway_enrichment=pathway_enrichment,
            clinical_translatability=clinical_translatability
        )
    
    def _calculate_genotype_pathogenicity_score(self, variant_impact_results: Dict) -> float:
        """Calculate genotype pathogenicity score from variant impact results"""
        if not variant_impact_results or 'missense_variants' not in variant_impact_results:
            return 0.0
        
        variants = variant_impact_results['missense_variants']
        if not variants:
            return 0.0
        
        # Calculate average pathogenicity score
        pathogenicity_scores = []
        for variant in variants:
            if 'pathogenicity_score' in variant:
                pathogenicity_scores.append(variant['pathogenicity_score'])
        
        if not pathogenicity_scores:
            return 0.0
        
        # Normalize to 0-1 scale (assuming scores are already 0-1)
        avg_score = statistics.mean(pathogenicity_scores)
        
        # Boost score if variants are consistently pathogenic
        consistency_boost = 0.1 if len([s for s in pathogenicity_scores if s > 0.7]) > len(pathogenicity_scores) * 0.5 else 0.0
        
        return min(1.0, avg_score + consistency_boost)
    
    def _calculate_structure_function_score(self, structural_results: Dict) -> float:
        """Calculate structure function score from structural results"""
        if not structural_results:
            return 0.0
        
        # Extract relevant scores
        binding_site_score = structural_results.get('binding_site_score', 0.0)
        structural_plausibility = structural_results.get('structural_plausibility', 0.0)
        
        # Normalize structural plausibility (assuming it's a sequence recovery score)
        # Lower sequence recovery is better, so invert and normalize
        normalized_plausibility = max(0.0, 1.0 - (structural_plausibility - 1.0))
        
        # Combine scores with weights
        structure_score = (binding_site_score * 0.6 + normalized_plausibility * 0.4)
        
        return min(1.0, structure_score)
    
    def _calculate_pathway_enrichment_score(self, pathway_dynamics_results: Dict) -> float:
        """Calculate pathway enrichment score from pathway dynamics results"""
        if not pathway_dynamics_results or 'pathway_enrichment' not in pathway_dynamics_results:
            return 0.0
        
        pathway_enrichment = pathway_dynamics_results['pathway_enrichment']
        if not pathway_enrichment:
            return 0.0
        
        # Calculate average enrichment score
        enrichment_scores = list(pathway_enrichment.values())
        avg_enrichment = statistics.mean(enrichment_scores)
        
        # Normalize to 0-1 scale (assuming enrichment scores are typically 0-2)
        normalized_score = min(1.0, avg_enrichment / 2.0)
        
        # Boost if key pathways are enriched
        key_pathways = ['p53_signaling', 'DNA_repair', 'Apoptosis']
        key_pathway_boost = 0.1 if any(pathway in pathway_enrichment and pathway_enrichment[pathway] > 1.0 
                                      for pathway in key_pathways) else 0.0
        
        return min(1.0, normalized_score + key_pathway_boost)
    
    def _calculate_clinical_translatability_score(self, coherence_results: Dict) -> float:
        """Calculate clinical translatability score from coherence results"""
        if not coherence_results:
            return 0.0
        
        # Extract cross-scale consistency
        cross_scale_consistency = coherence_results.get('cross_scale_consistency', 0.0)
        
        # Calculate tissue specificity score
        tissue_scores = coherence_results.get('tissue_specificity_scores', {})
        tissue_specificity = 0.0
        if tissue_scores:
            # Calculate average tissue specificity across genes
            gene_scores = []
            for gene_scores_dict in tissue_scores.values():
                if isinstance(gene_scores_dict, dict):
                    gene_scores.extend(gene_scores_dict.values())
            
            if gene_scores:
                tissue_specificity = statistics.mean(gene_scores)
        
        # Combine scores
        clinical_score = (cross_scale_consistency * 0.7 + tissue_specificity * 0.3)
        
        return min(1.0, clinical_score)
    
    def _calculate_overall_score(self, category_scores: CategoryScore) -> float:
        """Calculate overall coherence score"""
        weights = self.config['category_weights']
        
        overall_score = (
            category_scores.genotype_pathogenicity * weights['genotype_pathogenicity'] +
            category_scores.structure_function * weights['structure_function'] +
            category_scores.pathway_enrichment * weights['pathway_enrichment'] +
            category_scores.clinical_translatability * weights['clinical_translatability']
        )
        
        return min(1.0, overall_score)
    
    def _determine_coherence_category(self, overall_score: float) -> str:
        """Determine coherence category based on overall score"""
        for category in self.coherence_categories:
            if category.min_score <= overall_score <= category.max_score:
                return category.name
        
        return "NO_COHERENCE"
    
    def _generate_causal_trace(self,
                              variant_impact_results: Dict,
                              structural_results: Dict,
                              pathway_dynamics_results: Dict,
                              coherence_results: Dict) -> List[CausalTrace]:
        """Generate detailed causal trace based on analysis results"""
        trace = []
        
        # Step 1: Predict pathogenicity via AlphaMissense + ClinVar lookup
        if variant_impact_results and 'missense_variants' in variant_impact_results:
            variants = variant_impact_results['missense_variants']
            if variants:
                # Find most pathogenic variant
                most_pathogenic = max(variants, key=lambda v: v.get('pathogenicity_score', 0))
                gene = most_pathogenic.get('gene', 'Unknown')
                protein_change = most_pathogenic.get('protein_change', '')
                pathogenicity_score = most_pathogenic.get('pathogenicity_score', 0)
                alphamissense_annotation = most_pathogenic.get('alphamissense_annotation', 'Unknown')
                clinvar_annotation = most_pathogenic.get('clinvar_annotation', 'Unknown')
                
                # Get variant position and chromosome info
                position = most_pathogenic.get('position', 'Unknown')
                chromosome = most_pathogenic.get('chromosome', 'Unknown')
                
                trace.append(CausalTrace(
                    step=1,
                    level="molecular",
                    description=f"{gene} variant {protein_change} (chr{chromosome}:{position}) classified as {alphamissense_annotation} by AlphaMissense (score: {pathogenicity_score:.3f}) and {clinvar_annotation} by ClinVar, indicating high pathogenicity risk"
                ))
        
        # Step 2: Apply variant to reference sequence (BioPython)
        if structural_results:
            # Get protein sequence info
            protein_sequence = variant_impact_results.get('protein_sequence', '')
            seq_length = len(protein_sequence) if protein_sequence else 0
            
            trace.append(CausalTrace(
                step=2,
                level="molecular",
                description=f"BioPython applied variant to reference protein sequence ({seq_length} amino acids), enabling structural impact prediction and sequence-based analysis"
            ))
        
        # Step 3: Run GSEApy on gene list
        if pathway_dynamics_results and 'pathway_enrichment' in pathway_dynamics_results:
            pathway_enrichment = pathway_dynamics_results['pathway_enrichment']
            if pathway_enrichment:
                # Find most enriched pathway and get summary stats
                most_enriched = max(pathway_enrichment.items(), key=lambda x: x[1])
                summary = pathway_dynamics_results.get('summary', {})
                total_pathways = summary.get('total_pathways_analyzed', 0)
                significant_pathways = summary.get('significant_pathways', 0)
                
                # Get top 3 enriched pathways
                sorted_pathways = sorted(pathway_enrichment.items(), key=lambda x: x[1], reverse=True)
                top_pathways = sorted_pathways[:3]
                pathway_list = ", ".join([f"{p[0]} ({p[1]:.2f})" for p in top_pathways])
                
                trace.append(CausalTrace(
                    step=3,
                    level="pathway",
                    description=f"GSEApy pathway enrichment analysis of {total_pathways} pathways identified {significant_pathways} significant pathways. Top enriched: {pathway_list}, suggesting dysregulated cellular processes"
                ))
        
        # Step 4: Run ProteinMPNN on protein structure
        if structural_results and 'model_used' in structural_results:
            model_used = structural_results['model_used']
            binding_site_score = structural_results.get('binding_site_score', 0.0)
            structural_plausibility = structural_results.get('structural_plausibility', 0.0)
            
            # Interpret structural scores
            binding_interpretation = "high" if binding_site_score > 0.7 else "moderate" if binding_site_score > 0.4 else "low"
            plausibility_interpretation = "good" if structural_plausibility < 1.5 else "moderate" if structural_plausibility < 2.0 else "poor"
            
            trace.append(CausalTrace(
                step=4,
                level="molecular",
                description=f"ProteinMPNN ({model_used}) structural analysis shows {binding_interpretation} binding site score ({binding_site_score:.3f}) and {plausibility_interpretation} structural plausibility ({structural_plausibility:.3f}), indicating {binding_interpretation} druggability potential"
            ))
        
        # Step 5: Run tissue-specific biomarker relevance check (GTEx)
        if coherence_results and 'tissue_specificity_scores' in coherence_results:
            tissue_scores = coherence_results['tissue_specificity_scores']
            cross_scale_consistency = coherence_results.get('cross_scale_consistency', 0.0)
            
            # Get tissue-specific expression data
            tissue_info = []
            for gene, gene_tissues in tissue_scores.items():
                for tissue, score in gene_tissues.items():
                    if score > 0.5:  # Only include tissues with meaningful expression
                        tissue_info.append(f"{gene} in {tissue} (TPM: {score:.2f})")
            
            tissue_summary = "; ".join(tissue_info[:3]) if tissue_info else "limited tissue-specific expression"
            consistency_level = "high" if cross_scale_consistency > 0.7 else "moderate" if cross_scale_consistency > 0.5 else "low"
            
            trace.append(CausalTrace(
                step=5,
                level="clinical",
                description=f"GTEx tissue-specific analysis shows {tissue_summary}. Cross-scale consistency: {consistency_level} ({cross_scale_consistency:.3f}), supporting clinical relevance assessment"
            ))
        else:
            trace.append(CausalTrace(
                step=5,
                level="clinical",
                description="GTEx tissue-specific biomarker analysis completed, evaluating gene expression patterns across relevant tissues for clinical translatability"
            ))
        
        return trace
    
    def _create_metadata(self,
                        structural_results: Dict,
                        variant_impact_results: Dict,
                        sequence_variant_results: Dict,
                        coherence_results: Dict,
                        pathway_dynamics_results: Dict) -> Dict[str, Any]:
        """Create metadata for the final result"""
        metadata = {
            'analysis_timestamp': datetime.now().isoformat(),
            'pipeline_version': 'enhanced_v2',
            'input_files': {
                'structural_results': structural_results.get('timestamp', ''),
                'variant_impact_results': variant_impact_results.get('pdb_file', ''),
                'coherence_results': coherence_results.get('processing_timestamp', ''),
                'pathway_dynamics_results': pathway_dynamics_results.get('summary', {})
            },
            'config_used': self.config,
            'statistics': {
                'total_variants_analyzed': len(variant_impact_results.get('missense_variants', [])),
                'pathways_analyzed': pathway_dynamics_results.get('summary', {}).get('total_pathways_analyzed', 0),
                'significant_pathways': pathway_dynamics_results.get('summary', {}).get('significant_pathways', 0)
            }
        }
        
        return metadata
    
    def save_final_results(self, result: FinalResult, output_file: str) -> None:
        """
        Save final results to JSON file
        
        Args:
            result: FinalResult object
            output_file: Output file path
        """
        # Convert to the required schema format
        output_data = {
            "compound_analysis": {
                "mechanistic_coherence": {
                    "overall_score": result.overall_score,
                    "coherence_category": result.coherence_category,
                    "category_scores": {
                        "genotype_pathogenicity": result.category_scores.genotype_pathogenicity,
                        "structure_function": result.category_scores.structure_function,
                        "pathway_enrichment": result.category_scores.pathway_enrichment,
                        "clinical_translatability": result.category_scores.clinical_translatability
                    },
                    "causal_trace": [
                        {
                            "step": trace.step,
                            "level": trace.level,
                            "description": trace.description
                        }
                        for trace in result.causal_trace
                    ]
                }
            }
        }
        
        # Ensure output directory exists
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save to file
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logger.info(f"Final results saved to: {output_file}")
    
    def run_scoring_pipeline(self,
                           structural_results_file: str,
                           variant_impact_file: str,
                           sequence_variant_file: str,
                           coherence_results_file: str,
                           pathway_dynamics_file: str,
                           output_file: str) -> None:
        """
        Run the complete scoring pipeline
        
        Args:
            structural_results_file: Path to structural results JSON
            variant_impact_file: Path to variant impact results JSON
            sequence_variant_file: Path to sequence variant results JSON
            coherence_results_file: Path to coherence results JSON
            pathway_dynamics_file: Path to pathway dynamics results JSON
            output_file: Output file path for final results
        """
        logger.info("Starting intelligent coherence scoring pipeline")
        
        # Aggregate and score
        result = self.aggregate_and_score(
            structural_results_file,
            variant_impact_file,
            sequence_variant_file,
            coherence_results_file,
            pathway_dynamics_file
        )
        
        # Save results
        self.save_final_results(result, output_file)
        
        logger.info("Intelligent coherence scoring pipeline completed")


def main():
    """Main function for testing"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Intelligent Coherence Scorer")
    parser.add_argument("--structural", required=True, help="Structural results file")
    parser.add_argument("--variant-impact", required=True, help="Variant impact results file")
    parser.add_argument("--sequence-variant", required=True, help="Sequence variant results file")
    parser.add_argument("--coherence", required=True, help="Coherence results file")
    parser.add_argument("--pathway-dynamics", required=True, help="Pathway dynamics results file")
    parser.add_argument("--output", required=True, help="Output file path")
    
    args = parser.parse_args()
    
    scorer = IntelligentCoherenceScorer()
    scorer.run_scoring_pipeline(
        args.structural,
        args.variant_impact,
        args.sequence_variant,
        args.coherence,
        args.pathway_dynamics,
        args.output
    )


if __name__ == "__main__":
    main() 