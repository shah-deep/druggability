#!/usr/bin/env python3
"""
Druggability Scoring Engine for Binding Pocket Druggability Scorer
Task 3: Generate composite druggability scores for detected pockets

This script computes a composite score using feature descriptors, with an option to 
train a simple ML model based on benchmark data if provided.

Author: Joshua Robert
Date: July 2025
"""

import os
import sys
import json
import argparse
from typing import Dict, List
import numpy as np
import warnings
warnings.filterwarnings('ignore')

try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import LinearRegression
    SKLEARN_AVAILABLE = True
except ImportError:
    print("Warning: scikit-learn not available. Install with: pip install scikit-learn")
    SKLEARN_AVAILABLE = False


class ScoringEngine:
    """
    A class for computing druggability scores.
    """

    def __init__(self):
        # Default weights for a weighted sum model
        self.default_weights = {
            'volume': 0.2,
            'surface_area': 0.1,
            'hydrophobicity': 0.25,
            'polarity': -0.15,   # Negative weight (lower polarity preferred)
            'shape_complexity': 0.1,
            'solvent_accessibility': 0.1,
            'charge': -0.1       # Slight penalty for high net charge
        }

        self.scaler = None
        self.model = None

    def load_features(self, feature_json_file: str) -> Dict[str, Dict]:
        """Load pocket features from JSON."""
        try:
            with open(feature_json_file, 'r') as f:
                features = json.load(f)
            return features
        except Exception as e:
            raise RuntimeError(f"Error loading features: {e}")

    def compute_weighted_score(self, pocket_features: Dict) -> float:
        """Compute composite score using weighted sum."""
        score = 0.0
        for feature, weight in self.default_weights.items():
            value = pocket_features.get(feature, 0.0)
            score += weight * value
        return score

    def prepare_ml_data(self, features: Dict[str, Dict]) -> (np.ndarray, List[str]):
        """Convert feature dictionary to numpy array for ML model."""
        # Use selected numeric features
        selected_features = [
            'volume', 'surface_area', 'hydrophobicity', 'polarity', 
            'shape_complexity', 'solvent_accessibility', 'charge'
        ]

        data = []
        pocket_ids = []

        for p_id, feats in features.items():
            row = [feats.get(f, 0.0) for f in selected_features]
            data.append(row)
            pocket_ids.append(p_id)

        return np.array(data), pocket_ids

    def train_ml_model(self, X: np.ndarray, y: np.ndarray, model_type='rf'):
        """Train a simple ML model."""
        if not SKLEARN_AVAILABLE:
            raise RuntimeError("scikit-learn is required for ML model training")

        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(X)

        if model_type == 'rf':
            model = RandomForestRegressor(n_estimators=200, random_state=42)
        else:
            model = LinearRegression()

        model.fit(X_scaled, y)
        self.model = model

    def predict_ml_score(self, X: np.ndarray) -> np.ndarray:
        """Predict scores using trained ML model."""
        if self.model is None or self.scaler is None:
            raise RuntimeError("ML model has not been trained")

        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)

    def generate_scores(self, feature_json_file: str, output_json: str, 
                       benchmark_file: str = None, model_type='rf') -> None:
        """
        Main function to generate druggability scores.

        Args:
            feature_json_file (str): Feature JSON path
            output_json (str): Output JSON path
            benchmark_file (str): Optional benchmark CSV with known scores
            model_type (str): 'rf' or 'linear'
        """
        # Load features
        all_features = self.load_features(feature_json_file)

        if not all_features:
            raise RuntimeError("No features found in input file")

        # Prepare data for ML model
        X, pocket_ids = self.prepare_ml_data(all_features)

        # If benchmark provided, train ML model
        if benchmark_file and SKLEARN_AVAILABLE:
            print("Training ML model using benchmark data ...")
            # Expect CSV: pocket_id, score
            bench_data = {}
            with open(benchmark_file, 'r') as f:
                for line in f:
                    parts = line.strip().split(',')
                    if len(parts) == 2:
                        bench_data[parts[0]] = float(parts[1])

            # Match benchmark pockets
            y = []
            X_train = []
            for idx, p_id in enumerate(pocket_ids):
                if p_id in bench_data:
                    y.append(bench_data[p_id])
                    X_train.append(X[idx])

            if len(y) >= 10:  # Need at least 10 data points to train
                y = np.array(y)
                X_train = np.array(X_train)
                self.train_ml_model(X_train, y, model_type=model_type)
                ml_scores = self.predict_ml_score(X)
            else:
                print("Insufficient benchmark data, falling back to weighted sum model")
                ml_scores = None
        else:
            ml_scores = None

        # Compute scores
        output_data = {}
        for idx, p_id in enumerate(pocket_ids):
            feats = all_features[p_id]
            wsum_score = self.compute_weighted_score(feats)

            if ml_scores is not None:
                final_score = (wsum_score + ml_scores[idx]) / 2.0
            else:
                final_score = wsum_score

            output_data[p_id] = {
                'composite_score': round(final_score, 3),
                'weighted_sum_score': round(wsum_score, 3),
                'ml_score': round(ml_scores[idx], 3) if ml_scores is not None else None,
                'features': feats
            }

        # Sort pockets by score (descending)
        sorted_pockets = sorted(output_data.items(), 
                               key=lambda item: item[1]['composite_score'], 
                               reverse=True)
        output_data_sorted = {p_id: data for p_id, data in sorted_pockets}

        # Save to JSON
        with open(output_json, 'w') as f:
            json.dump(output_data_sorted, f, indent=2)

        print(f"Scores saved to {output_json}")
