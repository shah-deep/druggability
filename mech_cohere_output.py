"""
Mechanistic Coherence Output System
Formats results and integrates with existing pipeline, including both Task 1 (druggability) and Task 2 (mechanistic coherence) outputs.
"""

import numpy as np
import json
from typing import Dict, List
import logging
from mech_cohere_scorer import MechanisticCoherenceScorer
from typing import Optional

class MechanisticCoherenceOutput:
    """Handles output formatting and integration"""

    def __init__(self, scorer: MechanisticCoherenceScorer):
        self.scorer = scorer

    def generate_coherence_report(self, model_outputs: Dict, molecule_id: str = "unknown") -> Dict:
        """Generate comprehensive mechanistic coherence report"""

        # Calculate category scores
        category_scores = self.scorer.calculate_category_scores(model_outputs)

        # Calculate final score
        final_score = self.scorer.calculate_final_score(category_scores)

        # Detect conflicts
        conflicts = []
        for cat in model_outputs.values():
            if isinstance(cat, dict) and "flags" in cat:
                conflicts.extend(cat["flags"])
            elif isinstance(cat, list):
                for item in cat:
                    if isinstance(item, dict) and "flags" in item:
                        conflicts.extend(item["flags"])

        # Calculate confidence
        confidence = self.scorer.calculate_confidence(model_outputs, conflicts)

        # Generate summary
        summary = self._generate_mechanism_summary(category_scores, conflicts, model_outputs)

        # Expand model_outputs with detailed formatting
        model_outputs = self._format_model_outputs(model_outputs)

        # Create structured output
        coherence_report = {
            "molecule_id": molecule_id,
            "mechanistic_coherence_score": round(final_score, 3),
            "category_scores": {k: round(v, 3) for k, v in category_scores.items()},
            "confidence": round(confidence, 3),
            "conflict_flags": conflicts,
            "summary": summary,
            "model_outputs": model_outputs,
            "analysis_metadata": {
                "timestamp": self._get_timestamp(),
                "model_versions": self._get_model_versions(),
                "scoring_weights": self.scorer.weights
            }
        }

        return coherence_report

    def _format_model_outputs(self, model_outputs: Dict) -> Dict:
        """Ensure detailed outputs for all used models are present in a consistent format."""
        formatted = {}
        for category, models in model_outputs.items():
            if isinstance(models, list):
                detailed = []
                for m in models:
                    if isinstance(m, dict):
                        entry = {k: v for k, v in m.items() if k != "flags"}  # strip internal-only fields
                        detailed.append(entry)
                formatted[category] = detailed
            else:
                formatted[category] = models
        return formatted

    def _generate_mechanism_summary(self, category_scores: Dict, conflicts: List[str], model_outputs: Dict) -> str:
        """Generate mechanism-to-English summary"""
        parts = []

        if conflicts:
            parts.append(f"⚠️ Conflicts detected: {', '.join(conflicts)}.")
        else:
            parts.append("No model conflicts detected.")

        score = category_scores.get("structure", 0)
        if score >= 0.75:
            parts.append("Structural models (ProteinMPNN, RFDiffusion) suggest stable protein folding and compatible binding geometry.")
        elif score >= 0.5:
            parts.append("Structural analysis moderately supports binding and stability.")
        else:
            parts.append("Structure models show limited evidence for functional binding.")

        score = category_scores.get("genotype", 0)
        if score >= 0.75:
            parts.append("Genotype predictors indicate strong support for gene expression or splicing correction.")
        elif score >= 0.5:
            parts.append("Genomic predictors provide partial support for the mechanism.")
        else:
            parts.append("Genotype predictors show weak or inconsistent effects.")

        score = category_scores.get("simulation", 0)
        if score >= 0.75:
            parts.append("Simulation models (e.g., scFoundation) predict favorable downstream response.")
        elif score >= 0.5:
            parts.append("Simulation results suggest modest or cell-specific responses.")
        else:
            parts.append("Simulation models indicate minimal impact.")

        score = category_scores.get("biomarker", 0)
        if score >= 0.75:
            parts.append("Biomarker models validate clinical relevance or observability.")
        elif score >= 0.5:
            parts.append("Biomarker results provide limited clinical insight.")
        else:
            parts.append("No strong biomarker evidence detected.")

        score = category_scores.get("metabolic", 0)
        if score >= 0.75:
            parts.append("Metabolic pathway models suggest downstream functional viability.")
        elif score >= 0.5:
            parts.append("Metabolic effects are partially supported.")
        else:
            parts.append("Metabolic simulations show no clear downstream coherence.")

        return " ".join(parts)

    def _get_timestamp(self) -> str:
        from datetime import datetime
        return datetime.now().isoformat()

    def _get_model_versions(self) -> Dict:
        return {
            "enformer": "v1.0",
            "bigrna": "v1.0",
            "alphamissense": "v1.0",
            "proteinmpnn": "v1.0",
            "rfdiffusion": "v1.0",
            "cobrapy": "v0.27.0"
        }

    def integrate_with_druggability(
        self,
        coherence_report: Dict,
        druggability_data: Dict,
        pocket_features: Dict = None,
        druggability_scores: Dict = None,
        druggability_summary: str = None
    ) -> Dict:

        task1_block = {}

        # Always load from pocket_features_protein_42.json
        try:
            with open("pocket_features_protein_42.json") as f:
                pocket_features = json.load(f)
            task1_block["pockets"] = pocket_features

            # Generate synthetic scores from fpocket_druggability_score
            synthetic_scores = {}
            for pocket_id, features in pocket_features.items():
                score = features.get("fpocket_druggability_score", 0.0)
                synthetic_scores[pocket_id] = {
                    "features": features,
                    "composite_score": score,
                    "weighted_sum_score": score,
                    "ml_score": None
                }

            task1_block["scores"] = synthetic_scores
            if not druggability_summary:
                top = max(synthetic_scores.items(), key=lambda x: x[1]["composite_score"])
                task1_block["summary"] = f"Top scoring pocket is {top[0]} with a composite score of {round(top[1]['composite_score'], 3)}."
        except Exception as e:
            print("Warning: Failed to load or process pocket_features_protein_42.json:", e)
            task1_block["summary"] = "Druggability features could not be loaded."

        integrated_report = {
            "compound_analysis": {
                "druggability_assessment": task1_block,
                "mechanistic_coherence": coherence_report
            },
            "combined_score": self._calculate_combined_score(coherence_report, task1_block),
            "recommendation": self._generate_recommendation(coherence_report, task1_block),
            "analysis_timestamp": self._get_timestamp()
        }

        return integrated_report

    def _calculate_combined_score(self, coherence_report: Dict, druggability_data: Dict) -> Dict:
        coherence_score = coherence_report.get('mechanistic_coherence_score', 0.5)
        druggability_scores = []

        def extract_score(pocket: Dict) -> Optional[float]:
            # Check for flat structure
            if "fpocket_druggability_score" in pocket:
                return pocket["fpocket_druggability_score"]
            # Check nested "features" dictionary
            if "features" in pocket and isinstance(pocket["features"], dict):
                return pocket["features"].get("fpocket_druggability_score")
            return None

        if isinstance(druggability_data, dict):
            # Task 1 style structure: { "scores": { "Pocket_1": {...}, ... } }
            if "scores" in druggability_data:
                for pocket in druggability_data["scores"].values():
                    score = extract_score(pocket)
                    if score is not None:
                        druggability_scores.append(score)
            else:
                # Possibly directly passing pocket_features dict
                for pocket in druggability_data.values():
                    if isinstance(pocket, dict):
                        score = extract_score(pocket)
                        if score is not None:
                            druggability_scores.append(score)

        # Use max score found or fallback
        best_druggability = max(druggability_scores) if druggability_scores else 0.0
        combined_score = 0.6 * coherence_score + 0.4 * best_druggability

        return {
            "combined_score": round(combined_score, 3),
            "mechanistic_weight": 0.6,
            "druggability_weight": 0.4,
            "mechanistic_score": round(coherence_score, 3),
            "best_druggability_score": round(best_druggability, 3)
        }


    def _generate_recommendation(self, coherence_report: Dict, druggability_data: Dict) -> str:
        coherence_score = coherence_report.get('mechanistic_coherence_score', 0.5)
        conflicts = len(coherence_report.get('conflict_flags', []))
        confidence = coherence_report.get('confidence', 0.5)

        if coherence_score > 0.7 and confidence > 0.8 and conflicts == 0:
            return "Highly recommended for further development"
        elif coherence_score > 0.6 and confidence > 0.7:
            return "Recommended with additional validation"
        elif coherence_score > 0.4:
            return "Consider with caution - requires mechanism optimization"
        else:
            return "Not recommended - poor mechanistic coherence"
