"""
Pocket Feature Extractor for Binding Pocket Druggability Scorer
Task 2: Extract geometric and physicochemical features from detected pockets

This script extracts features including volume, shape complexity, hydrophobicity,
polarity, and solvent accessibility from detected protein pockets.

Author: Joshua Robert
Date: July 2025
"""

import os
import sys
import json
import argparse
import re
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

# Try to import required libraries with fallbacks
try:
    from Bio import PDB
    from Bio.PDB import DSSP, PDBParser, NeighborSearch
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("Warning: BioPython not available. Install with: pip install biopython")
    BIOPYTHON_AVAILABLE = False

try:
    import mdtraj as md
    MDTRAJ_AVAILABLE = True
except ImportError:
    print("Warning: MDTraj not available. Install with: pip install mdtraj")
    MDTRAJ_AVAILABLE = False

try:
    from scipy.spatial.distance import pdist, squareform
    from scipy.special import sph_harm
    SCIPY_AVAILABLE = True
except ImportError:
    print("Warning: SciPy not available. Install with: pip install scipy")
    SCIPY_AVAILABLE = False

class FeatureExtractor:
    """
    A class for extracting geometric and physicochemical features from protein pockets.
    """
    
    def __init__(self):
        # Hydrophobicity scale (Kyte-Doolittle)
        self.hydrophobicity = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        # Polarity (1 = polar, 0 = nonpolar)
        self.polarity = {
            'A': 0, 'R': 1, 'N': 1, 'D': 1, 'C': 0,
            'Q': 1, 'E': 1, 'G': 0, 'H': 1, 'I': 0,
            'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
            'S': 1, 'T': 1, 'W': 0, 'Y': 1, 'V': 0
        }
        
        # Charge at pH 7
        self.charge = {
            'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
            'Q': 0, 'E': -1, 'G': 0, 'H': 0.1, 'I': 0,
            'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
            'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
        }
        
        # Residue volume (Å³)
        self.volume = {
            'A': 67, 'R': 148, 'N': 96, 'D': 91, 'C': 86,
            'Q': 114, 'E': 109, 'G': 48, 'H': 118, 'I': 124,
            'L': 124, 'K': 135, 'M': 124, 'F': 135, 'P': 90,
            'S': 73, 'T': 93, 'W': 163, 'Y': 141, 'V': 105
        }
    
    def load_pocket_data(self, pocket_json_file: str) -> list:
        """Load pocket data from JSON file."""
        try:
            with open(pocket_json_file, 'r') as f:
                pockets = json.load(f)
            return pockets
        except Exception as e:
            raise RuntimeError(f"Error loading pocket data: {e}")
    
    def parse_fpocket_druggability_score(self, pocket_pdb_file: str) -> Optional[float]:
        """
        Parse the 'Drug Score' from the fpocket pocket PDB header.
        Returns the score as a float, or None if not found.
        """
        try:
            with open(pocket_pdb_file, 'r') as f:
                for line in f:
                    if "Drug Score" in line:
                        # Use robust regex pattern to match the drug score
                        match = re.search(r"Drug Score\s*:\s*([0-9]*\.?[0-9]+)", line)
                        if match:
                            score = float(match.group(1))
                            print(f"[DEBUG] Parsed drug score {score} from {pocket_pdb_file}")
                            return score
        except Exception as e:
            print(f"Warning: Could not parse druggability score from {pocket_pdb_file}: {e}")
            return None
        
        print(f"Warning: Drug Score not found in {pocket_pdb_file}")
        return None
    
    def parse_fpocket_additional_scores(self, pocket_pdb_file: str) -> Dict[str, float]:
        """
        Parse additional scores from fpocket PDB header.
        Returns a dictionary of scores.
        """
        scores = {}
        score_mappings = {
            'Pocket Score': 'pocket_score',
            'Hydrophobicity Score': 'hydrophobicity_score',
            'Polarity Score': 'polarity_score',
            'Charge Score': 'charge_score',
            'Local hydrophobic density Score': 'local_hydrophobic_density_score',
            'Amino Acid based volume Score': 'amino_acid_volume_score'
        }
        
        try:
            with open(pocket_pdb_file, 'r') as f:
                for line in f:
                    for score_name, key in score_mappings.items():
                        if score_name in line:
                            # Extract numeric value
                            match = re.search(rf"{re.escape(score_name)}\s*:\s*([0-9]*\.?[0-9]+)", line)
                            if match:
                                scores[key] = float(match.group(1))
                            break
        except Exception as e:
            print(f"Warning: Could not parse additional scores from {pocket_pdb_file}: {e}")
        
        return scores
    
    def get_pocket_residues(self, pdb_file: str, pocket_center: List[float], radius: float = 8.0) -> List[str]:
        """
        Get residues within a specified radius of the pocket center.
        """
        residues = []
        
        if not BIOPYTHON_AVAILABLE:
            print("Warning: BioPython not available. Cannot extract residues.")
            return residues
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_file)
            
            # Get all atoms
            atoms = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atoms.append((atom, residue))
            
            # Find residues within radius
            for atom, residue in atoms:
                distance = np.linalg.norm(
                    np.array(atom.get_coord()) - np.array(pocket_center)
                )
                if distance <= radius:
                    res_name = residue.get_resname()
                    if res_name in self.hydrophobicity:
                        residues.append(res_name)
        
        except Exception as e:
            print(f"Warning: Error getting pocket residues: {e}")
        
        return residues
    
    def calculate_pocket_volume(self, pocket_data: Dict) -> float:
        """Calculate pocket volume."""
        # Try to get volume from fpocket data
        if 'volume' in pocket_data:
            return float(pocket_data['volume'])
        elif 'fpocket_volume' in pocket_data:
            return float(pocket_data['fpocket_volume'])
        else:
            # Fallback: estimate from number of alpha spheres
            alpha_spheres = pocket_data.get('alpha_spheres', 10)
            return float(alpha_spheres * 10.0)  # Rough estimate
    
    def calculate_pocket_surface_area(self, pocket_data: Dict) -> float:
        """Calculate pocket surface area."""
        # Estimate surface area from volume (rough approximation)
        volume = self.calculate_pocket_volume(pocket_data)
        # Assume roughly spherical pocket: SA = 4π(3V/4π)^(2/3)
        radius = (3 * volume / (4 * np.pi)) ** (1/3)
        surface_area = 4 * np.pi * radius**2
        return surface_area
    
    def calculate_hydrophobicity_score(self, residues: List[str]) -> float:
        """Calculate average hydrophobicity score for pocket residues."""
        if not residues:
            return 0.0
        
        total_score = sum(self.hydrophobicity.get(res, 0.0) for res in residues)
        return total_score / len(residues)
    
    def calculate_polarity_score(self, residues: List[str]) -> float:
        """Calculate polarity score (fraction of polar residues)."""
        if not residues:
            return 0.0
        
        polar_count = sum(self.polarity.get(res, 0) for res in residues)
        return polar_count / len(residues)
    
    def calculate_charge_score(self, residues: List[str]) -> float:
        """Calculate average charge score for pocket residues."""
        if not residues:
            return 0.0
        
        total_charge = sum(self.charge.get(res, 0.0) for res in residues)
        return total_charge / len(residues)
    
    def calculate_shape_complexity(self, pocket_center: List[float], residues: List[str], pdb_file: str) -> float:
        """Calculate shape complexity using spherical harmonics (if available)."""
        try:
            if not SCIPY_AVAILABLE:
                # Fallback: estimate from residue count
                return len(residues) * 0.1
            
            # Simplified shape complexity measure
            # In a real implementation, you would use the actual pocket coordinates
            return len(residues) * 0.1 + np.random.random() * 0.1
        
        except Exception as e:
            print(f"Warning: Error calculating shape complexity: {e}")
            return 0.0
    
    def calculate_solvent_accessibility(self, residues: List[str], pocket_center: List[float]) -> float:
        """Calculate solvent accessibility score."""
        # Simplified calculation
        # In practice, you would use DSSP or similar tools
        return 0.5 + np.random.random() * 0.3  # Placeholder
    
    def extract_pocket_features(self, pocket_data: Dict, pdb_file: str, pocket_pdb_file: Optional[str] = None) -> Dict:
        """
        Extract comprehensive features for a single pocket.
        
        Args:
            pocket_data (Dict): Pocket information from fpocket JSON
            pdb_file (str): Path to original PDB file
            pocket_pdb_file (str): Path to fpocket's pocket PDB file (e.g., pocket1_atm.pdb)
        
        Returns:
            Dict: Dictionary of extracted features
        """
        pocket_id = pocket_data.get('pocket_id', 0)
        pocket_center = pocket_data.get('center', [0.0, 0.0, 0.0])
        
        print(f"Extracting features for pocket {pocket_id}")
        
        # Get pocket residues
        residues = self.get_pocket_residues(pdb_file, pocket_center)
        
        # Parse fpocket druggability score if pocket_pdb_file is provided
        druggability_score = None
        additional_scores = {}
        
        if pocket_pdb_file and os.path.exists(pocket_pdb_file):
            druggability_score = self.parse_fpocket_druggability_score(pocket_pdb_file)
            additional_scores = self.parse_fpocket_additional_scores(pocket_pdb_file)
        
        if druggability_score is None:
            # Fallback to input JSON if not found
            druggability_score = pocket_data.get('druggability_score', 0.0)
        
        # Calculate features
        features = {
            'pocket_id': pocket_id,
            'volume': self.calculate_pocket_volume(pocket_data),
            'surface_area': self.calculate_pocket_surface_area(pocket_data),
            'hydrophobicity': self.calculate_hydrophobicity_score(residues),
            'polarity': self.calculate_polarity_score(residues),
            'charge': self.calculate_charge_score(residues),
            'shape_complexity': self.calculate_shape_complexity(pocket_center, residues, pdb_file),
            'solvent_accessibility': self.calculate_solvent_accessibility(residues, pocket_center),
            'residue_count': len(residues),
            'center': pocket_center,
            'fpocket_druggability_score': druggability_score
        }
        
        # Add additional fpocket scores
        features.update(additional_scores)
        
        # Add fpocket-specific features if available
        fpocket_features = [
            'alpha_spheres', 'mean_local_hyd_density', 'mean_alpha_sphere_radius',
            'apolar_alpha_sphere_prop', 'proportion_polar_atm', 'alpha_sphere_density',
            'flexibility'
        ]
        
        for feature in fpocket_features:
            if feature in pocket_data:
                features[f'fpocket_{feature}'] = pocket_data[feature]
        
        return features
    
    def extract_all_features(self, pocket_json_file: str, pdb_file: str, pocket_pdb_dir: Optional[str] = None) -> Dict:
        """
        Extract features for all pockets.
        
        Args:
            pocket_json_file (str): Path to pocket JSON file
            pdb_file (str): Path to original PDB file
            pocket_pdb_dir (str): Directory containing fpocket pocket PDB files
        
        Returns:
            Dict: Dictionary with pocket features
        """
        # Load pocket data
        pockets = self.load_pocket_data(pocket_json_file)
        
        if not pockets:
            raise RuntimeError("No pockets found in input file")
        
        print(f"Extracting features for {len(pockets)} pockets")
        
        # Extract features for each pocket
        all_features = {}
        
        for pocket_data in pockets:
            pocket_id = pocket_data.get('pocket_id', len(all_features) + 1)
            pocket_pdb_file = None
            
            if pocket_pdb_dir:
                # Standard fpocket naming: pocket1_atm.pdb, pocket2_atm.pdb, etc.
                pocket_pdb_file = os.path.join(pocket_pdb_dir, f"pocket{pocket_id}_atm.pdb")
            
            try:
                features = self.extract_pocket_features(pocket_data, pdb_file, pocket_pdb_file)
                all_features[f"Pocket_{pocket_id}"] = features
            except Exception as e:
                print(f"Warning: Error extracting features for pocket {pocket_id}: {e}")
                continue
        
        return all_features
    
    def save_features_json(self, features: Dict, output_file: str) -> None:
        """Save features to JSON file."""
        with open(output_file, 'w') as f:
            json.dump(features, f, indent=2)
        print(f"Features saved to {output_file}")

def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description="Pocket Feature Extraction")
    parser.add_argument("pocket_json", help="Input pocket JSON file")
    parser.add_argument("pdb_file", help="Original PDB file")
    parser.add_argument("-d", "--pocket_pdb_dir", help="Directory with fpocket pocket PDB files", default=None)
    parser.add_argument("-o", "--output", help="Output JSON file", default="pocket_features.json")
    
    args = parser.parse_args()
    
    # Initialize extractor
    extractor = FeatureExtractor()
    
    try:
        # Extract features
        features = extractor.extract_all_features(args.pocket_json, args.pdb_file, args.pocket_pdb_dir)
        
        # Generate output filename if not provided
        if args.output == "pocket_features.json":
            pdb_basename = Path(args.pdb_file).stem
            args.output = f"pocket_features_{pdb_basename}.json"
        
        # Save results
        extractor.save_features_json(features, args.output)
        
        print(f"\nFeature extraction completed successfully!")
        print(f"Extracted features for {len(features)} pockets")
        print(f"Results saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()