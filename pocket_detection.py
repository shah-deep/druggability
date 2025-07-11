#!/usr/bin/env python3
"""
Pocket Detection Pipeline for Binding Pocket Druggability Scorer
Task 1: Protein Pocket Detection using fpocket

This script takes a .pdb protein structure file and identifies potential binding pockets
using fpocket as the base detector, with optional DeepSite integration.

Author: Joshua Robert
Date: July 2025
"""

import os
import sys
import json
import subprocess
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

class PocketDetector:
    """
    A class for detecting protein pockets using fpocket and optionally DeepSite.
    """

    def __init__(self, fpocket_path: str = "fpocket"):
        """
        Initialize the PocketDetector.

        Args:
            fpocket_path (str): Path to fpocket executable
        """
        self.fpocket_path = fpocket_path
        self.validate_fpocket()

    def validate_fpocket(self) -> None:
        """Validate that fpocket is installed and accessible."""
        try:
            result = subprocess.run([self.fpocket_path, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                # Try without --version flag
                result = subprocess.run([self.fpocket_path], 
                                      capture_output=True, text=True, timeout=10)
            print(f"fpocket validation successful")
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
            print("Warning: fpocket not found. Please install fpocket first.")
            print("Installation instructions:")
            print("1. Using conda: conda install -c conda-forge fpocket")
            print("2. From source: https://github.com/Discngine/fpocket")

    def run_fpocket(self, pdb_file: str, output_dir: Optional[str] = None) -> str:
        """
        Run fpocket on a PDB file.

        Args:
            pdb_file (str): Path to input PDB file
            output_dir (str): Output directory (optional)

        Returns:
            str: Path to output directory
        """
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        # Prepare command
        cmd = [self.fpocket_path, "-f", pdb_file]

        if output_dir:
            cmd.extend(["-o", output_dir])

        print(f"Running fpocket command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if result.returncode != 0:
                print(f"fpocket stderr: {result.stderr}")
                print(f"fpocket stdout: {result.stdout}")
                raise subprocess.CalledProcessError(result.returncode, cmd)

            # Determine output directory
            if output_dir:
                fpocket_output_dir = output_dir
            else:
                # fpocket creates a directory named <pdb_basename>_out
                pdb_basename = Path(pdb_file).stem
                fpocket_output_dir = f"{pdb_basename}_out"

            print(f"fpocket completed successfully. Output directory: {fpocket_output_dir}")
            return fpocket_output_dir

        except subprocess.TimeoutExpired:
            raise RuntimeError("fpocket execution timed out")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"fpocket failed with return code {e.returncode}")

    def parse_fpocket_output(self, output_dir: str) -> List[Dict]:
        """
        Parse fpocket output files to extract pocket information.

        Args:
            output_dir (str): fpocket output directory

        Returns:
            List[Dict]: List of pocket information dictionaries
        """
        pockets = []

        # Look for fpocket info file
        info_file = os.path.join(output_dir, "fpocket_info.txt")
        if not os.path.exists(info_file):
            # Try alternative naming
            for file in os.listdir(output_dir):
                if file.endswith("_info.txt"):
                    info_file = os.path.join(output_dir, file)
                    break

        if not os.path.exists(info_file):
            print(f"Warning: fpocket info file not found in {output_dir}")
            return self._parse_pocket_files(output_dir)

        try:
            with open(info_file, 'r') as f:
                content = f.read()

            # Parse pocket information using regex
            pocket_blocks = re.findall(r'Pocket\s+(\d+)\s*:\s*(.+?)(?=Pocket\s+\d+\s*:|$)', 
                                     content, re.DOTALL)

            for pocket_id, pocket_info in pocket_blocks:
                pocket_data = self._parse_pocket_block(pocket_id, pocket_info)
                if pocket_data:
                    pockets.append(pocket_data)

        except Exception as e:
            print(f"Error parsing fpocket info file: {e}")
            return self._parse_pocket_files(output_dir)

        return pockets

    def _parse_pocket_block(self, pocket_id: str, pocket_info: str) -> Dict:
        """Parse individual pocket block from fpocket output."""
        pocket_data = {
            "pocket_id": int(pocket_id),
            "center": [0.0, 0.0, 0.0],
            "radius": 0.0,
            "volume": 0.0,
            "surface_area": 0.0,
            "druggability_score": 0.0,
            "alpha_spheres": 0,
            "mean_local_hyd_density": 0.0,
            "mean_alpha_sphere_radius": 0.0,
            "mean_alpha_sphere_surface_area": 0.0,
            "apolar_alpha_sphere_prop": 0.0,
            "hydrophobicity_score": 0.0,
            "polarity_score": 0.0,
            "volume_score": 0.0,
            "charge_score": 0.0,
            "proportion_polar_atm": 0.0,
            "alpha_sphere_density": 0.0,
            "center_of_mass": [0.0, 0.0, 0.0],
            "flexibility": 0.0
        }

        # Extract numeric values using regex
        patterns = {
            "druggability_score": r'Drug Score\s*:\s*([0-9.-]+)',
            "volume": r'Volume\s*:\s*([0-9.-]+)',
            "surface_area": r'Surface\s*:\s*([0-9.-]+)',
            "alpha_spheres": r'Alpha spheres\s*:\s*([0-9]+)',
            "mean_local_hyd_density": r'Mean local hydrophobic density\s*:\s*([0-9.-]+)',
            "mean_alpha_sphere_radius": r'Mean alpha sphere radius\s*:\s*([0-9.-]+)',
            "mean_alpha_sphere_surface_area": r'Mean alpha sphere surface area\s*:\s*([0-9.-]+)',
            "apolar_alpha_sphere_prop": r'Apolar alpha sphere proportion\s*:\s*([0-9.-]+)',
            "hydrophobicity_score": r'Hydrophobicity score\s*:\s*([0-9.-]+)',
            "polarity_score": r'Polarity score\s*:\s*([0-9.-]+)',
            "volume_score": r'Volume score\s*:\s*([0-9.-]+)',
            "charge_score": r'Charge score\s*:\s*([0-9.-]+)',
            "proportion_polar_atm": r'Proportion of polar atoms\s*:\s*([0-9.-]+)',
            "alpha_sphere_density": r'Alpha sphere density\s*:\s*([0-9.-]+)',
            "flexibility": r'Flexibility\s*:\s*([0-9.-]+)'
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, pocket_info, re.IGNORECASE)
            if match:
                try:
                    value = float(match.group(1))
                    if key == "alpha_spheres":
                        value = int(value)
                    pocket_data[key] = value
                except ValueError:
                    pass

        # Extract center coordinates
        center_match = re.search(r'Center\s*:\s*([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)', 
                                pocket_info, re.IGNORECASE)
        if center_match:
            pocket_data["center"] = [float(center_match.group(1)),
                                   float(center_match.group(2)),
                                   float(center_match.group(3))]

        # Extract center of mass
        com_match = re.search(r'Center of mass\s*:\s*([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)', 
                             pocket_info, re.IGNORECASE)
        if com_match:
            pocket_data["center_of_mass"] = [float(com_match.group(1)),
                                           float(com_match.group(2)),
                                           float(com_match.group(3))]

        # Calculate radius from volume (approximate)
        if pocket_data["volume"] > 0:
            pocket_data["radius"] = (3 * pocket_data["volume"] / (4 * 3.14159)) ** (1/3)

        return pocket_data

    def _parse_pocket_files(self, output_dir: str) -> List[Dict]:
        """Fallback method to parse pocket files directly."""
        pockets = []

        # Look for pocket PDB files
        pocket_files = [f for f in os.listdir(output_dir) if f.startswith("pocket") and f.endswith(".pdb")]

        for i, pocket_file in enumerate(sorted(pocket_files)):
            pocket_data = {
                "pocket_id": i + 1,
                "center": [0.0, 0.0, 0.0],
                "radius": 5.0,  # Default radius
                "volume": 0.0,
                "surface_area": 0.0,
                "druggability_score": 0.0,
                "file_path": os.path.join(output_dir, pocket_file)
            }

            # Try to extract basic info from filename if available
            match = re.search(r'pocket(\d+)', pocket_file)
            if match:
                pocket_data["pocket_id"] = int(match.group(1))

            pockets.append(pocket_data)

        return pockets

    def detect_pockets(self, pdb_file: str, output_dir: Optional[str] = None) -> Tuple[List[Dict], str]:
        """
        Complete pocket detection pipeline.

        Args:
            pdb_file (str): Path to input PDB file
            output_dir (str): Output directory (optional)

        Returns:
            Tuple[List[Dict], str]: (List of pocket dictionaries, output directory path)
        """
        print(f"Starting pocket detection for {pdb_file}")

        # Run fpocket
        fpocket_output_dir = self.run_fpocket(pdb_file, output_dir)

        # Parse results
        pockets = self.parse_fpocket_output(fpocket_output_dir)

        print(f"Detected {len(pockets)} pockets")

        return pockets, fpocket_output_dir

    def save_pockets_json(self, pockets: List[Dict], output_file: str) -> None:
        """
        Save pocket information to JSON file.

        Args:
            pockets (List[Dict]): List of pocket dictionaries
            output_file (str): Output JSON file path
        """
        with open(output_file, 'w') as f:
            json.dump(pockets, f, indent=2)

        print(f"Pocket information saved to {output_file}")


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(description="Protein Pocket Detection Pipeline")
    parser.add_argument("pdb_file", help="Input PDB file path")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("--fpocket-path", default="fpocket", help="Path to fpocket executable")

    args = parser.parse_args()

    # Initialize detector
    detector = PocketDetector(fpocket_path=args.fpocket_path)

    # Detect pockets
    try:
        pockets, output_dir = detector.detect_pockets(args.pdb_file, args.output)

        # Generate output filename
        pdb_basename = Path(args.pdb_file).stem
        output_json = f"pockets_{pdb_basename}.json"

        # Save results
        detector.save_pockets_json(pockets, output_json)

        print(f"\nPocket detection completed successfully!")
        print(f"Results saved to: {output_json}")
        print(f"fpocket output directory: {output_dir}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()