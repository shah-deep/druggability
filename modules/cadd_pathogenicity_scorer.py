#!/usr/bin/env python3
"""
CADD Pathogenicity Scorer - Gold Standard Variant Pathogenicity Assessment

This module implements CADD (Combined Annotation Dependent Depletion) scoring,
which is the most widely accepted and standardized method for variant pathogenicity assessment.

CADD integrates multiple annotations including:
- Conservation scores
- Regulatory annotations  
- Protein function predictions
- Disease associations
- Population frequencies

Author: Enhanced for ITS4.3
Date: January 2025
"""

import os
import re
import json
import logging
import requests
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, field
from datetime import datetime
import time
import sqlite3
import gzip
import hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.parse
from gene_coordinate_mapper import GeneCoordinateMapper

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class CADDScore:
    """CADD pathogenicity score result"""
    variant_id: str
    cadd_raw_score: float
    cadd_phred_score: float
    chromosome: str
    position: int
    reference_allele: str
    alternate_allele: str
    gene: Optional[str] = None
    transcript: Optional[str] = None
    protein_change: Optional[str] = None
    consequence: Optional[str] = None
    confidence: float = 1.0
    prediction_source: str = "CADD"
    warnings: List[str] = field(default_factory=list)
    
    @property
    def pathogenicity_level(self) -> str:
        """Determine pathogenicity level based on CADD Phred score"""
        if self.cadd_phred_score >= 30:
            return "HIGH"
        elif self.cadd_phred_score >= 20:
            return "MEDIUM"
        elif self.cadd_phred_score >= 10:
            return "LOW"
        else:
            return "BENIGN"
    
    @property
    def is_pathogenic(self) -> bool:
        """Determine if variant is likely pathogenic (CADD Phred >= 20)"""
        return self.cadd_phred_score >= 20

class CADDClient:
    """Client for CADD pathogenicity scoring via public CADD API (bihealth.org)"""
    
    def __init__(self, cache_dir: str = "cache/cadd", use_api: bool = True):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.api_version = "GRCh38-v1.7"  # Use latest public version
        self.coordinate_mapper = GeneCoordinateMapper()
        logger.info(f"CADD Client initialized with public API: {use_api}")

    def score_variants(self, variants: List[Dict]) -> Dict[str, CADDScore]:
        """Score variants using CADD public API with coordinate mapping"""
        results = {}
        for variant in variants:
            variant_id = variant.get('id', 'unknown')
            gene = variant.get('gene', '')
            protein_change = variant.get('protein_change', '')
            chrom = variant.get('chromosome', variant.get('chr', ''))
            pos = variant.get('position', variant.get('pos', 0))
            ref = variant.get('reference_allele', variant.get('ref', ''))
            alt = variant.get('alternate_allele', variant.get('alt', ''))
            
            # Map to CADD coordinates if needed
            if not (chrom and pos and ref and alt):
                if gene and protein_change:
                    logger.info(f"Mapping {gene} {protein_change} to genomic coordinates")
                    coordinates = self.coordinate_mapper.map_variant_to_coordinates(gene, protein_change, chrom, pos)
                    chrom = coordinates['chromosome']
                    pos = coordinates['position']
                    ref = coordinates['reference_allele']
                    alt = coordinates['alternate_allele']
                    logger.info(f"Mapped to {chrom}:{pos}_{ref}_{alt}")
                else:
                    raise ValueError(f"Insufficient information for variant {variant_id}")
            
            # Only support SNVs
            if not (len(ref) == 1 and len(alt) == 1):
                raise ValueError(f"Only SNVs are supported by the CADD API. Variant {variant_id} is not an SNV.")
            
            url = f"https://cadd.gs.washington.edu/api/v1.0/{self.api_version}/{chrom}:{pos}_{ref}_{alt}"
            logger.info(f"Querying CADD API for {variant_id}: {url}")
            resp = requests.get(url, timeout=30)
            if resp.status_code != 200:
                raise RuntimeError(f"CADD API error for {variant_id}: {resp.status_code} - {resp.text}")
            data = resp.json()
            if not data or not isinstance(data, list) or len(data) == 0:
                raise ValueError(f"No CADD score found for {variant_id} ({chrom}:{pos}_{ref}_{alt})")
            entry = data[0]
            cadd_score = CADDScore(
                variant_id=variant_id,
                cadd_raw_score=float(entry['RawScore']),
                cadd_phred_score=float(entry['PHRED']),
                chromosome=chrom,
                position=int(pos),
                reference_allele=ref,
                alternate_allele=alt,
                gene=variant.get('gene'),
                protein_change=variant.get('protein_change'),
                prediction_source="CADD Public API"
            )
            results[variant_id] = cadd_score
        return results
    
    def _init_cache_db(self):
        """Initialize SQLite cache database"""
        db_path = self.cache_dir / "cadd_cache.db"
        self.db_path = db_path
        
        with sqlite3.connect(db_path) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS cadd_scores (
                    variant_id TEXT PRIMARY KEY,
                    cadd_raw_score REAL,
                    cadd_phred_score REAL,
                    chromosome TEXT,
                    position INTEGER,
                    reference_allele TEXT,
                    alternate_allele TEXT,
                    gene TEXT,
                    transcript TEXT,
                    protein_change TEXT,
                    consequence TEXT,
                    timestamp TEXT
                )
            """)
            conn.commit()
    

    
    def _get_cached_scores(self, variants: List[Dict]) -> Dict[str, CADDScore]:
        """Get cached CADD scores from database"""
        cached_results = {}
        
        with sqlite3.connect(self.db_path) as conn:
            for variant in variants:
                variant_id = variant.get('id', 'unknown')
                
                # Create cache key
                cache_key = self._create_cache_key(variant)
                
                cursor = conn.execute("""
                    SELECT cadd_raw_score, cadd_phred_score, chromosome, position,
                           reference_allele, alternate_allele, gene, transcript,
                           protein_change, consequence, timestamp
                    FROM cadd_scores WHERE variant_id = ?
                """, (cache_key,))
                
                row = cursor.fetchone()
                if row:
                    cadd_score = CADDScore(
                        variant_id=variant_id,
                        cadd_raw_score=row[0],
                        cadd_phred_score=row[1],
                        chromosome=row[2],
                        position=row[3],
                        reference_allele=row[4],
                        alternate_allele=row[5],
                        gene=row[6],
                        transcript=row[7],
                        protein_change=row[8],
                        consequence=row[9],
                        prediction_source="CADD (cached)"
                    )
                    cached_results[variant_id] = cadd_score
                    
        return cached_results
    
    def _cache_scores(self, scores: Dict[str, CADDScore]):
        """Cache CADD scores in database"""
        with sqlite3.connect(self.db_path) as conn:
            for variant_id, score in scores.items():
                cache_key = self._create_cache_key({
                    'chromosome': score.chromosome,
                    'position': score.position,
                    'reference_allele': score.reference_allele,
                    'alternate_allele': score.alternate_allele
                })
                
                conn.execute("""
                    INSERT OR REPLACE INTO cadd_scores 
                    (variant_id, cadd_raw_score, cadd_phred_score, chromosome, position,
                     reference_allele, alternate_allele, gene, transcript, protein_change,
                     consequence, timestamp)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    cache_key, score.cadd_raw_score, score.cadd_phred_score,
                    score.chromosome, score.position, score.reference_allele,
                    score.alternate_allele, score.gene, score.transcript,
                    score.protein_change, score.consequence, datetime.now().isoformat()
                ))
            conn.commit()
    
    def _create_cache_key(self, variant: Dict) -> str:
        """Create unique cache key for variant"""
        # Use chromosome:position:ref:alt format
        chrom = variant.get('chromosome', variant.get('chr', ''))
        pos = variant.get('position', variant.get('pos', ''))
        ref = variant.get('reference_allele', variant.get('ref', ''))
        alt = variant.get('alternate_allele', variant.get('alt', ''))
        
        return f"{chrom}:{pos}:{ref}:{alt}"
    
    def _query_cadd_api(self, variants: List[Dict]) -> Dict[str, CADDScore]:
        """Query CADD API for pathogenicity scores"""
        results = {}
        
        # Prepare variants for API
        api_variants = []
        for variant in variants:
            # Parse variant information
            chrom = variant.get('chromosome', variant.get('chr', ''))
            pos = variant.get('position', variant.get('pos', ''))
            ref = variant.get('reference_allele', variant.get('ref', ''))
            alt = variant.get('alternate_allele', variant.get('alt', ''))
            
            if not all([chrom, pos, ref, alt]):
                logger.warning(f"Skipping variant with incomplete information: {variant}")
                continue
                
            api_variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt
            })
        
        if not api_variants:
            return results
        
        try:
            # Use the correct CADD API endpoint
            api_url = "https://cadd.gs.washington.edu/api/v1.0/annotations"
            
            # Prepare the request payload in the correct format
            # CADD API expects variants in the format: chr:pos:ref:alt
            payload = {
                'variants': [
                    f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}"
                    for v in api_variants
                ]
            }
            
            # Try the main API endpoint
            try:
                response = requests.post(
                    api_url,
                    json=payload,
                    headers={'Content-Type': 'application/json'},
                    timeout=30
                )
            except:
                # Try alternative endpoint
                api_url = "https://cadd.gs.washington.edu/api/v1.0/snv"
                response = requests.post(
                    api_url,
                    json=payload,
                    headers={'Content-Type': 'application/json'},
                    timeout=30
                )
            
            if response.status_code == 200:
                data = response.json()
                
                for i, variant in enumerate(variants):
                    variant_id = variant.get('id', f'variant_{i}')
                    
                    if i < len(data.get('variants', [])):
                        variant_data = data['variants'][i]
                        
                        # Extract CADD scores
                        cadd_raw = variant_data.get('cadd_raw', 0.0)
                        cadd_phred = variant_data.get('cadd_phred', 0.0)
                        
                        # Extract additional information
                        chrom = variant_data.get('chrom', '')
                        pos = variant_data.get('pos', 0)
                        ref = variant_data.get('ref', '')
                        alt = variant_data.get('alt', '')
                        
                        # Extract gene and transcript information
                        gene = None
                        transcript = None
                        protein_change = None
                        consequence = None
                        
                        if 'annotations' in variant_data and variant_data['annotations']:
                            annotations = variant_data['annotations']
                            if annotations and len(annotations) > 0:
                                # Get first annotation (usually most relevant)
                                ann = annotations[0]
                                if ann:
                                    gene = ann.get('Gene_Name')
                                    transcript = ann.get('Feature_ID')
                                    protein_change = ann.get('Protein_Change')
                                    consequence = ann.get('Consequence')
                        
                        cadd_score = CADDScore(
                            variant_id=variant_id,
                            cadd_raw_score=cadd_raw,
                            cadd_phred_score=cadd_phred,
                            chromosome=chrom,
                            position=pos,
                            reference_allele=ref,
                            alternate_allele=alt,
                            gene=gene,
                            transcript=transcript,
                            protein_change=protein_change,
                            consequence=consequence,
                            prediction_source="CADD API"
                        )
                        
                        results[variant_id] = cadd_score
                    else:
                        logger.warning(f"No CADD data returned for variant {variant_id}")
                        
            else:
                logger.error(f"CADD API error: {response.status_code} - {response.text}")
                raise RuntimeError(f"CADD API failed with status {response.status_code}: {response.text}")
                
        except Exception as e:
            logger.error(f"Error querying CADD API: {e}")
            raise RuntimeError(f"Failed to query CADD API: {e}")
            
        return results
    


class PathogenicityScorer:
    """Main class for variant pathogenicity scoring using CADD"""
    
    def __init__(self, config: Optional[Dict] = None):
        self.config = config or {}
        
        # Initialize CADD client
        use_api = self.config.get('use_cadd_api', True)
        cache_dir = self.config.get('cadd_cache_dir', 'cache/cadd')
        
        self.cadd_client = CADDClient(cache_dir=cache_dir, use_api=use_api)
        
        logger.info("Pathogenicity Scorer initialized with CADD integration")
    
    def score_variants(self, variants: List[Dict]) -> Dict[str, CADDScore]:
        """Score variants for pathogenicity using CADD"""
        logger.info(f"Scoring {len(variants)} variants for pathogenicity")
        
        # Validate input variants
        validated_variants = self._validate_variants(variants)
        
        if not validated_variants:
            raise ValueError("No valid variants provided for scoring")
        
        # Score variants using CADD
        results = self.cadd_client.score_variants(validated_variants)
        
        logger.info(f"Successfully scored {len(results)} variants")
        return results
    
    def _validate_variants(self, variants: List[Dict]) -> List[Dict]:
        """Validate and clean variant input"""
        validated = []
        
        for variant in variants:
            # Check required fields
            required_fields = ['chromosome', 'position', 'reference_allele', 'alternate_allele']
            missing_fields = [field for field in required_fields if not variant.get(field)]
            
            if missing_fields:
                logger.warning(f"Skipping variant with missing fields {missing_fields}: {variant}")
                continue
            
            # Normalize chromosome format
            chrom = variant.get('chromosome', variant.get('chr', ''))
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            variant['chromosome'] = chrom
            
            # Ensure position is integer
            try:
                variant['position'] = int(variant.get('position', variant.get('pos', 0)))
            except (ValueError, TypeError):
                logger.warning(f"Invalid position for variant: {variant}")
                continue
            
            validated.append(variant)
        
        return validated
    
    def get_summary_statistics(self, results: Dict[str, CADDScore]) -> Dict[str, Any]:
        """Generate summary statistics for scored variants"""
        if not results:
            return {
                'total_variants': 0,
                'mean_cadd_phred': 0.0,
                'median_cadd_phred': 0.0,
                'std_cadd_phred': 0.0,
                'min_cadd_phred': 0.0,
                'max_cadd_phred': 0.0,
                'pathogenicity_distribution': {
                    'HIGH': 0,
                    'MEDIUM': 0,
                    'LOW': 0,
                    'BENIGN': 0
                },
                'pathogenic_variants': 0
            }
        
        scores = [score.cadd_phred_score for score in results.values()]
        
        summary = {
            'total_variants': len(results),
            'mean_cadd_phred': np.mean(scores) if scores else 0.0,
            'median_cadd_phred': np.median(scores) if scores else 0.0,
            'std_cadd_phred': np.std(scores) if scores else 0.0,
            'min_cadd_phred': np.min(scores) if scores else 0.0,
            'max_cadd_phred': np.max(scores) if scores else 0.0,
            'pathogenicity_distribution': {
                'HIGH': len([s for s in results.values() if s.pathogenicity_level == 'HIGH']),
                'MEDIUM': len([s for s in results.values() if s.pathogenicity_level == 'MEDIUM']),
                'LOW': len([s for s in results.values() if s.pathogenicity_level == 'LOW']),
                'BENIGN': len([s for s in results.values() if s.pathogenicity_level == 'BENIGN'])
            },
            'pathogenic_variants': len([s for s in results.values() if s.is_pathogenic])
        }
        
        return summary
    
    def export_results(self, results: Dict[str, CADDScore], output_path: str):
        """Export CADD scoring results to file"""
        output_data = {
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'total_variants': len(results),
                'scoring_method': 'CADD',
                'summary': self.get_summary_statistics(results)
            },
            'variants': {}
        }
        
        for variant_id, score in results.items():
            output_data['variants'][variant_id] = {
                'cadd_raw_score': score.cadd_raw_score,
                'cadd_phred_score': score.cadd_phred_score,
                'pathogenicity_level': score.pathogenicity_level,
                'is_pathogenic': score.is_pathogenic,
                'chromosome': score.chromosome,
                'position': score.position,
                'reference_allele': score.reference_allele,
                'alternate_allele': score.alternate_allele,
                'gene': score.gene,
                'transcript': score.transcript,
                'protein_change': score.protein_change,
                'consequence': score.consequence,
                'prediction_source': score.prediction_source,
                'warnings': score.warnings
            }
        
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logger.info(f"Results exported to {output_path}")

def main():
    """Example usage of CADD Pathogenicity Scorer"""
    
    # Example variants
    variants = [
        {
            'id': 'rs121913343',
            'chromosome': '7',
            'position': 140453136,
            'reference_allele': 'A',
            'alternate_allele': 'T',
            'gene': 'BRAF'
        },
        {
            'id': 'rs121913342', 
            'chromosome': '7',
            'position': 140453135,
            'reference_allele': 'T',
            'alternate_allele': 'A',
            'gene': 'BRAF'
        }
    ]
    
    # Initialize scorer
    scorer = PathogenicityScorer()
    
    # Score variants
    results = scorer.score_variants(variants)
    
    # Print results
    for variant_id, score in results.items():
        print(f"\nVariant: {variant_id}")
        print(f"CADD Phred Score: {score.cadd_phred_score:.2f}")
        print(f"Pathogenicity Level: {score.pathogenicity_level}")
        print(f"Is Pathogenic: {score.is_pathogenic}")
        print(f"Gene: {score.gene}")
        print(f"Protein Change: {score.protein_change}")
    
    # Export results
    scorer.export_results(results, 'cadd_pathogenicity_results.json')
    
    # Print summary
    summary = scorer.get_summary_statistics(results)
    print(f"\nSummary Statistics:")
    print(f"Total Variants: {summary['total_variants']}")
    print(f"Mean CADD Phred: {summary['mean_cadd_phred']:.2f}")
    print(f"Pathogenic Variants: {summary['pathogenic_variants']}")

if __name__ == "__main__":
    main() 