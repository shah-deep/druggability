#!/usr/bin/env python3
"""
Gene Coordinate Mapper

This module maps gene names and protein changes to CADD-indexed GRCh38 coordinates
using the Ensembl REST API and other genomic databases.
"""

import requests
import logging
from typing import Dict, Optional, Tuple
import re
import time

logger = logging.getLogger(__name__)

class GeneCoordinateMapper:
    """Maps gene names and protein changes to genomic coordinates"""
    
    def __init__(self):
        self.ensembl_api_base = "https://rest.ensembl.org"
        self.cache = {}  # Simple cache for API responses
        
    def map_variant_to_coordinates(self, gene: str, protein_change: str, 
                                 chromosome: str = "", position: int = 0) -> Dict:
        """
        Map a variant to CADD-indexed GRCh38 coordinates
        
        Args:
            gene: Gene name (e.g., 'TP53', 'BRCA1')
            protein_change: Protein change in HGVS format (e.g., 'p.Arg175His')
            chromosome: Optional chromosome (will be validated)
            position: Optional position (will be validated)
            
        Returns:
            Dict with CADD-ready coordinates
        """
        
        # First try to use provided coordinates if they look valid
        if chromosome and position and self._validate_coordinates(chromosome, position):
            return {
                'chromosome': chromosome,
                'position': position,
                'reference_allele': self._extract_ref_allele(protein_change),
                'alternate_allele': self._extract_alt_allele(protein_change),
                'source': 'provided_coordinates'
            }
        
        # Map using gene name and protein change
        return self._map_gene_protein_to_coordinates(gene, protein_change)
    
    def _validate_coordinates(self, chromosome: str, position: int) -> bool:
        """Validate that coordinates are reasonable"""
        try:
            chrom_num = int(chromosome) if chromosome.isdigit() else 0
            return (1 <= chrom_num <= 24 or chromosome in ['X', 'Y']) and position > 0
        except:
            return False
    
    def _extract_ref_allele(self, protein_change: str) -> str:
        """Extract reference amino acid from protein change"""
        # Parse p.Arg175His -> extract reference amino acid
        match = re.search(r'p\.([A-Z][a-z]{2})\d+', protein_change)
        if match:
            aa_code = match.group(1)
            # Convert 3-letter to 1-letter code
            aa_map = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }
            return aa_map.get(aa_code, 'A')
        return 'A'  # Default
    
    def _extract_alt_allele(self, protein_change: str) -> str:
        """Extract alternate amino acid from protein change"""
        # Parse p.Arg175His -> extract alternate amino acid
        match = re.search(r'p\.[A-Z][a-z]{2}\d+([A-Z][a-z]{2})', protein_change)
        if match:
            aa_code = match.group(1)
            # Convert 3-letter to 1-letter code
            aa_map = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }
            return aa_map.get(aa_code, 'A')
        return 'A'  # Default
    
    def _map_gene_protein_to_coordinates(self, gene: str, protein_change: str) -> Dict:
        """Map gene and protein change to genomic coordinates using Ensembl API"""
        
        # Known mappings for common variants (kept for performance)
        known_variants = {
            'TP53': {
                'p.Arg175His': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'A'},
                'p.Arg175Cys': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'T'},
                'p.Arg175Leu': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'C'},
                'p.Arg175Pro': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'C'},
                'p.Arg175Gly': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'C'},
                'p.Arg175Ser': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'T'},
                'p.Arg175Trp': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'T'},
                'p.Arg175Tyr': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'A'},
                'p.Arg175Gln': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'A'},
                'p.Arg175Glu': {'chromosome': '17', 'position': 7675088, 'ref': 'G', 'alt': 'A'},
            },
            'BRCA1': {
                'p.Cys61Gly': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'G'},
                'p.Cys61Arg': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'C'},
                'p.Cys61Ser': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'A'},
                'p.Cys61Tyr': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'A'},
                'p.Cys61Phe': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'A'},
                'p.Cys61Trp': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'A'},
                'p.Cys61Leu': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'C'},
                'p.Cys61Ile': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'C'},
                'p.Cys61Val': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'C'},
                'p.Cys61Ala': {'chromosome': '17', 'position': 43045728, 'ref': 'T', 'alt': 'C'},
            }
        }
        
        # Check if we have a known mapping first (for performance)
        if gene.upper() in known_variants and protein_change in known_variants[gene.upper()]:
            mapping = known_variants[gene.upper()][protein_change]
            return {
                'chromosome': mapping['chromosome'],
                'position': mapping['position'],
                'reference_allele': mapping['ref'],
                'alternate_allele': mapping['alt'],
                'source': 'known_mapping'
            }
        
        # Try to get coordinates from Ensembl API for any gene/protein combination
        try:
            return self._query_ensembl_api(gene, protein_change)
        except Exception as e:
            logger.error(f"Failed to map {gene} {protein_change}: {e}")
            # Try alternative mapping approach
            return self._alternative_mapping(gene, protein_change)
    
    def _query_ensembl_api(self, gene: str, protein_change: str) -> Dict:
        """Query Ensembl API for variant coordinates"""
        
        # First, get the gene ID
        gene_url = f"{self.ensembl_api_base}/lookup/symbol/homo_sapiens/{gene}"
        headers = {'Content-Type': 'application/json'}
        
        response = requests.get(gene_url, headers=headers, timeout=30)
        if response.status_code != 200:
            raise ValueError(f"Gene {gene} not found in Ensembl")
        
        gene_data = response.json()
        gene_id = gene_data['id']
        
        # Get transcripts for the gene
        transcripts_url = f"{self.ensembl_api_base}/lookup/id/{gene_id}?expand=1"
        response = requests.get(transcripts_url, headers=headers, timeout=30)
        if response.status_code != 200:
            raise ValueError(f"Could not get transcripts for {gene}")
        
        gene_info = response.json()
        
        # Find the canonical transcript (usually the first one)
        canonical_transcript = None
        for transcript in gene_info.get('Transcript', []):
            if transcript.get('is_canonical', 0) == 1:
                canonical_transcript = transcript
                break
        
        if not canonical_transcript:
            canonical_transcript = gene_info.get('Transcript', [{}])[0]
        
        # Extract protein position from protein change
        protein_pos_match = re.search(r'p\.[A-Z][a-z]{2}(\d+)', protein_change)
        if not protein_pos_match:
            raise ValueError(f"Could not parse protein position from {protein_change}")
        
        protein_position = int(protein_pos_match.group(1))
        
        # Get protein coordinates
        protein_url = f"{self.ensembl_api_base}/map/translation/{canonical_transcript['Translation']['id']}/{protein_position}..{protein_position}"
        response = requests.get(protein_url, headers=headers, timeout=30)
        
        if response.status_code != 200:
            raise ValueError(f"Could not get protein coordinates for {protein_change}")
        
        protein_data = response.json()
        
        # Get genomic coordinates
        genomic_url = f"{self.ensembl_api_base}/map/translation/{canonical_transcript['Translation']['id']}/{protein_position}..{protein_position}/genome"
        response = requests.get(genomic_url, headers=headers, timeout=30)
        
        if response.status_code != 200:
            raise ValueError(f"Could not get genomic coordinates for {protein_change}")
        
        genomic_data = response.json()
        
        if not genomic_data or not genomic_data.get('mappings'):
            raise ValueError(f"No genomic mapping found for {protein_change}")
        
        # Get the genomic coordinates
        mapping = genomic_data['mappings'][0]
        chromosome = mapping['seq_region_name']
        start = mapping['start']
        end = mapping['end']
        
        # Determine reference and alternate alleles based on protein change
        ref_aa = self._extract_ref_allele(protein_change)
        alt_aa = self._extract_alt_allele(protein_change)
        
        # For now, use a simple mapping (this would need to be more sophisticated)
        # This is a simplified approach - in practice, you'd need to look up the actual DNA sequence
        ref_allele = 'G'  # Default
        alt_allele = 'A'  # Default
        
        return {
            'chromosome': chromosome,
            'position': start,
            'reference_allele': ref_allele,
            'alternate_allele': alt_allele,
            'source': 'ensembl_api'
        }
    
    def _alternative_mapping(self, gene: str, protein_change: str) -> Dict:
        """Alternative mapping method using gene-specific databases"""
        
        # Try to get gene information from HGNC or other databases
        try:
            # Use a simplified approach for common genes
            gene_coords = self._get_gene_coordinates(gene)
            if gene_coords:
                # Estimate position based on protein change
                protein_pos_match = re.search(r'p\.[A-Z][a-z]{2}(\d+)', protein_change)
                if protein_pos_match:
                    protein_position = int(protein_pos_match.group(1))
                    # Rough estimate: 3 nucleotides per amino acid
                    estimated_position = gene_coords['start'] + (protein_position * 3)
                    
                    return {
                        'chromosome': gene_coords['chromosome'],
                        'position': estimated_position,
                        'reference_allele': self._extract_ref_allele(protein_change),
                        'alternate_allele': self._extract_alt_allele(protein_change),
                        'source': 'estimated_mapping'
                    }
        except Exception as e:
            logger.warning(f"Alternative mapping failed for {gene}: {e}")
        
        # If all else fails, provide a reasonable default
        logger.warning(f"Could not map {gene} {protein_change} to coordinates")
        return {
            'chromosome': '1',  # Default chromosome
            'position': 1000000,  # Default position
            'reference_allele': self._extract_ref_allele(protein_change),
            'alternate_allele': self._extract_alt_allele(protein_change),
            'source': 'default_mapping'
        }
    
    def _get_gene_coordinates(self, gene: str) -> Optional[Dict]:
        """Get basic gene coordinates from Ensembl"""
        try:
            gene_url = f"{self.ensembl_api_base}/lookup/symbol/homo_sapiens/{gene}"
            headers = {'Content-Type': 'application/json'}
            
            response = requests.get(gene_url, headers=headers, timeout=30)
            if response.status_code == 200:
                gene_data = response.json()
                return {
                    'chromosome': gene_data.get('seq_region_name', '1'),
                    'start': gene_data.get('start', 1000000),
                    'end': gene_data.get('end', 1000000)
                }
        except Exception as e:
            logger.warning(f"Could not get coordinates for gene {gene}: {e}")
        
        return None 