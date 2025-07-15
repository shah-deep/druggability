#!/usr/bin/env python3
"""
Variant Impact Analyzer - Parallel Version

A module for analyzing variant impact using AlphaMissense and ClinVar databases.
Performs parallel processing of variants using separate processes,
with asynchronous I/O within each process.

Example usage:
    analyzer = VariantImpactAnalyzer()
    results = analyzer.analyze_variants("processed_input.json")
"""

import os
import json
import sqlite3
import logging
import re
import asyncio
import multiprocessing as mp
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import statistics
from collections import Counter
import argparse
import pickle
import time
import hashlib

# Import the ClinVar annotator
try:
    from .clinvar_annotator import ClinVarAnnotator, ClinVarAnnotation
except ImportError:
    from clinvar_annotator import ClinVarAnnotator, ClinVarAnnotation

# Configure logging
from .logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)


@dataclass
class AlphaMissenseResult:
    """AlphaMissense annotation result"""
    variant_id: str
    gene: str
    protein_change: str
    average_pathogenicity_score: Optional[float] = None
    max_occurring_class: Optional[str] = None
    matching_transcripts: List[str] = field(default_factory=list)
    total_matches: int = 0
    confidence: str = "low"  # "high", "medium", "low", "none"
    warnings: List[str] = field(default_factory=list)


@dataclass
class VariantImpactResult:
    """Complete variant impact analysis result"""
    variant_id: str
    gene: str
    protein_change: str
    alphamissense: AlphaMissenseResult
    clinvar: ClinVarAnnotation
    processing_timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class AlphaMissenseCache:
    """Cache for AlphaMissense results to avoid repeated database queries"""
    
    def __init__(self, cache_dir: str = "cache"):
        """
        Initialize AlphaMissense cache
        
        Args:
            cache_dir: Directory to store cache database
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if the cache_dir already contains 'alphamissense' to avoid duplication
        if 'alphamissense' in str(self.cache_dir):
            self.cache_db_path = self.cache_dir / "alphamissense_cache.db"
        else:
            self.cache_db_path = self.cache_dir / "alphamissense" / "alphamissense_cache.db"
            self.cache_db_path.parent.mkdir(parents=True, exist_ok=True)
        
        self._init_cache_db()

    def _init_cache_db(self):
        """Initialize SQLite cache database"""
        try:
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS alphamissense_cache (
                    cache_key TEXT PRIMARY KEY,
                    position INTEGER,
                    reference TEXT,
                    alternate TEXT,
                    gene TEXT,
                    protein_change TEXT,
                    pathogenicity_score REAL,
                    alphamissense_annotation TEXT,
                    alphamissense_confidence TEXT,
                    matching_transcripts TEXT,
                    total_matches INTEGER,
                    warnings TEXT,
                    timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Create index for faster lookups
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_variant_lookup 
                ON alphamissense_cache(position, reference, alternate, gene, protein_change)
            """)
            
            conn.commit()
            conn.close()

            
        except Exception as e:
            logger.error(f"Error initializing AlphaMissense cache database: {e}")
    
    def _generate_cache_key(self, variant_data: Dict) -> str:
        """Generate a unique cache key for a variant"""
        # Create a hash of the key variant fields
        key_fields = [
            str(variant_data.get('position', '')),
            str(variant_data.get('reference', '')),
            str(variant_data.get('alternate', '')),
            str(variant_data.get('gene', '')),
            str(variant_data.get('protein_change', ''))
        ]
        key_string = '|'.join(key_fields)
        return hashlib.md5(key_string.encode()).hexdigest()
    
    def get_cached_result(self, variant_data: Dict) -> Optional[AlphaMissenseResult]:
        """
        Get cached AlphaMissense result for a variant
        
        Args:
            variant_data: Variant dictionary containing position, reference, alternate, gene, protein_change
            
        Returns:
            AlphaMissenseResult if found in cache, None otherwise
        """
        try:
            cache_key = self._generate_cache_key(variant_data)
            
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                SELECT pathogenicity_score, alphamissense_annotation, alphamissense_confidence,
                       matching_transcripts, total_matches, warnings
                FROM alphamissense_cache 
                WHERE cache_key = ?
            """, (cache_key,))
            
            result = cursor.fetchone()
            conn.close()
            
            if result:
                pathogenicity_score, annotation, confidence, transcripts_str, total_matches, warnings_str = result
                
                # Parse stored data
                matching_transcripts = transcripts_str.split('|') if transcripts_str else []
                warnings = warnings_str.split('|') if warnings_str else []
                return AlphaMissenseResult(
                    variant_id=variant_data.get('id', ''),
                    gene=variant_data.get('gene', ''),
                    protein_change=variant_data.get('protein_change', ''),
                    average_pathogenicity_score=pathogenicity_score,
                    max_occurring_class=annotation,
                    matching_transcripts=matching_transcripts,
                    total_matches=total_matches or 0,
                    confidence=confidence or "low",
                    warnings=warnings
                )
            
            return None
            
        except Exception as e:
            logger.error(f"Error retrieving from AlphaMissense cache: {e}")
            return None
    
    def cache_result(self, variant_data: Dict, result: AlphaMissenseResult):
        """
        Cache AlphaMissense result for a variant
        
        Args:
            variant_data: Variant dictionary
            result: AlphaMissenseResult to cache
        """
        try:
            cache_key = self._generate_cache_key(variant_data)
            
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            # Prepare data for storage
            matching_transcripts_str = '|'.join(result.matching_transcripts) if result.matching_transcripts else ''
            warnings_str = '|'.join(result.warnings) if result.warnings else ''
            
            cursor.execute("""
                INSERT OR REPLACE INTO alphamissense_cache 
                (cache_key, position, reference, alternate, gene, protein_change,
                 pathogenicity_score, alphamissense_annotation, alphamissense_confidence,
                 matching_transcripts, total_matches, warnings)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                cache_key,
                variant_data.get('position'),
                variant_data.get('reference'),
                variant_data.get('alternate'),
                variant_data.get('gene'),
                variant_data.get('protein_change'),
                result.average_pathogenicity_score,
                result.max_occurring_class,
                result.confidence,
                matching_transcripts_str,
                result.total_matches,
                warnings_str
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logger.error(f"Error caching AlphaMissense result: {e}")
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        try:
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute("SELECT COUNT(*) FROM alphamissense_cache")
            total_entries = cursor.fetchone()[0]
            
            cursor.execute("SELECT COUNT(DISTINCT gene) FROM alphamissense_cache")
            unique_genes = cursor.fetchone()[0]
            
            cursor.execute("SELECT AVG(pathogenicity_score) FROM alphamissense_cache WHERE pathogenicity_score IS NOT NULL")
            avg_score = cursor.fetchone()[0]
            
            # Count entries with "none" confidence (these are not cached)
            cursor.execute("SELECT COUNT(*) FROM alphamissense_cache WHERE alphamissense_confidence = 'none'")
            none_confidence_count = cursor.fetchone()[0]
            
            conn.close()
            
            return {
                'total_entries': total_entries,
                'unique_genes': unique_genes,
                'average_pathogenicity_score': avg_score,
                'cache_size_mb': self.cache_db_path.stat().st_size / (1024 * 1024) if self.cache_db_path.exists() else 0,
                'none_confidence_entries': none_confidence_count
            }
            
        except Exception as e:
            logger.error(f"Error getting cache stats: {e}")
            return {}


class VariantImpactAnalyzer:
    """Parallel analyzer for variant impact using AlphaMissense and ClinVar"""
    
    def __init__(self, alphamissense_db_path: str = "cache/alphamissense/alphamissense_hg38.db", 
                 max_workers: Optional[int] = None, cache_dir: str = "cache"):
        """
        Initialize the variant impact analyzer
        
        Args:
            alphamissense_db_path: Path to AlphaMissense database
            max_workers: Maximum number of parallel processes (default: CPU count)
            cache_dir: Directory for cache storage
        """
        self.alphamissense_db_path = Path(alphamissense_db_path)
        self.max_workers = max_workers or mp.cpu_count()
        self.cache = AlphaMissenseCache(cache_dir)
        
        # Check if database exists, if not, use cache-only mode
        if not self.alphamissense_db_path.exists():
            logger.warning(f"AlphaMissense database not found: {self.alphamissense_db_path}")
            logger.warning("Falling back to cache-only mode. Only previously cached variants will be processed.")
            self.cache_only_mode = True
        else:
            self.cache_only_mode = False
        
        logger.info(f"Variant Impact Analyzer initialized with {self.max_workers} workers")
        if self.cache_only_mode:
            logger.info("Running in cache-only mode")
    
    def analyze_variants(self, input_file: str) -> Dict[str, Any]:
        """
        Analyze all variants in the processed input file using parallel processing

        Args:
            input_file: Path to processed_input.json file

        Returns:
            Dictionary containing analysis results and original data
        """
        # Load input data
        with open(input_file, 'r') as f:
            raw_data = json.load(f)

        # Handle both list format and dictionary format
        if isinstance(raw_data, list):
            variants = raw_data
            input_data = {'missense_variants': variants}
        else:
            variants = raw_data.get('missense_variants', [])
            input_data = raw_data if isinstance(raw_data, dict) else {'missense_variants': variants}

        logger.info(f"Analyzing {len(variants)} variants using {self.max_workers} parallel processes")

        # Create temporary directory for inter-process communication
        temp_dir = Path(tempfile.mkdtemp(prefix="variant_analysis_"))

        try:
            # Process variants in parallel
            results = self._process_variants_parallel(variants, temp_dir)

            # Update variants with results
            for variant in variants:
                variant_id = variant['id']
                if variant_id in results:
                    result = results[variant_id]

                    # Add AlphaMissense results to variant
                    variant['pathogenicity_score'] = result.alphamissense.average_pathogenicity_score
                    variant['alphamissense_annotation'] = result.alphamissense.max_occurring_class.title() if result.alphamissense.max_occurring_class else None
                    variant['alphamissense_confidence'] = result.alphamissense.confidence
                    
                    # Add ClinVar results to variant
                    variant['clinvar_variation_id'] = result.clinvar.variation_id
                    # Default clinvar_annotation to "Uncertain_Significance" if not present or None/empty
                    clinvar_ann = result.clinvar.clinical_significance
                    if clinvar_ann is None or clinvar_ann == "":
                        variant['clinvar_annotation'] = "Uncertain_Significance"
                    else:
                        variant['clinvar_annotation'] = clinvar_ann.title()

            # Add analysis metadata
            # input_data['variant_impact_analysis'] = {  # type: ignore
            #     'analysis_timestamp': datetime.now().isoformat(),
            #     'total_variants': len(variants),
            #     'parallel_workers': self.max_workers,
            #     'cache_stats': self.cache.get_cache_stats()
            # }

            return input_data

        finally:
            # Clean up temporary directory
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    def _process_variants_parallel(self, variants: List[Dict], temp_dir: Path) -> Dict[str, VariantImpactResult]:
        """
        Process variants in parallel using multiprocessing
        
        Args:
            variants: List of variant dictionaries
            temp_dir: Temporary directory for inter-process communication
            
        Returns:
            Dictionary mapping variant_id to VariantImpactResult
        """
        # Prepare arguments for each process
        process_args = []
        for i, variant in enumerate(variants):
            # Create unique temp file for this variant
            temp_file = temp_dir / f"variant_{i}_{variant['id']}.pkl"
            process_args.append((variant, str(self.alphamissense_db_path), str(temp_file), str(self.cache.cache_db_path)))
        
        # Process variants in parallel
        with mp.Pool(processes=self.max_workers) as pool:
            pool.map(self._process_single_variant_worker, process_args)
        
        # Collect results from temp files
        results = {}
        for i, variant in enumerate(variants):
            temp_file = temp_dir / f"variant_{i}_{variant['id']}.pkl"
            if temp_file.exists():
                try:
                    with open(temp_file, 'rb') as f:
                        result = pickle.load(f)
                    results[variant['id']] = result
                except Exception as e:
                    logger.error(f"Error loading result for {variant['id']}: {e}")
        
        return results
    
    def save_results(self, results: Dict[str, Any], output_file: str):
        """
        Save analysis results to file
        
        Args:
            results: Analysis results dictionary
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        
    
    @staticmethod
    def _process_single_variant_worker(args: Tuple[Dict, str, str, str]) -> None:
        """
        Worker function to process a single variant in a separate process
        
        Args:
            args: Tuple of (variant, alphamissense_db_path, temp_file_path, cache_db_path)
        """
        variant, alphamissense_db_path, temp_file_path, cache_db_path = args
        variant_id = variant['id']
        
        try:
            # Create analyzer instance for this process
            analyzer = SingleVariantAnalyzer(alphamissense_db_path, cache_db_path)
            
            # Process variant asynchronously
            result = asyncio.run(analyzer.analyze_variant_async(variant))
            
            # Save result to temp file
            with open(temp_file_path, 'wb') as f:
                pickle.dump(result, f)
                

            
        except Exception as e:
            logger.error(f"Error processing variant {variant_id}: {e}")
            # Create error result
            error_result = VariantImpactResult(
                variant_id=variant_id,
                gene=variant.get('gene', ''),
                protein_change=variant.get('protein_change', ''),
                alphamissense=AlphaMissenseResult(
                    variant_id=variant_id,
                    gene=variant.get('gene', ''),
                    protein_change=variant.get('protein_change', ''),
                    warnings=[f"Processing error: {str(e)}"]
                ),
                clinvar=ClinVarAnnotation(
                    variant_id=variant_id,
                    gene=variant.get('gene', ''),
                    protein_change=variant.get('protein_change', ''),
                    clinical_significance="Uncertain_Significance",
                    warnings=[f"Processing error: {str(e)}"]
                )
            )
            with open(temp_file_path, 'wb') as f:
                pickle.dump(error_result, f)


class SingleVariantAnalyzer:
    """Analyzer for a single variant with async I/O"""
    
    def __init__(self, alphamissense_db_path: str, cache_db_path: str):
        """
        Initialize single variant analyzer
        
        Args:
            alphamissense_db_path: Path to AlphaMissense database
            cache_db_path: Path to cache database
        """
        self.alphamissense_db_path = Path(alphamissense_db_path)
        self.cache_db_path = Path(cache_db_path)
        self.clinvar_annotator = ClinVarAnnotator()
        # Use the cache directory (parent of the cache db file) to avoid path duplication
        cache_dir = self.cache_db_path.parent
        self.cache = AlphaMissenseCache(str(cache_dir))
        self.cache.cache_db_path = self.cache_db_path  # Override cache path for this process
        
        # Check if we're in cache-only mode
        self.cache_only_mode = not self.alphamissense_db_path.exists()
    
    async def analyze_variant_async(self, variant: Dict) -> VariantImpactResult:
        """
        Analyze a single variant using async I/O for AlphaMissense and ClinVar
        
        Args:
            variant: Variant dictionary from input
            
        Returns:
            VariantImpactResult with analysis results
        """
        variant_id = variant['id']
        gene = variant['gene']
        protein_change = variant['protein_change']
        vep_transcript_ids = variant.get('vep_transcript_ids', [])
        
        # Run AlphaMissense and ClinVar analysis concurrently
        alphamissense_task = asyncio.create_task(
            self._analyze_alphamissense_async(variant_id, gene, protein_change, vep_transcript_ids, variant)
        )
        clinvar_task = asyncio.create_task(
            self._analyze_clinvar_async(variant)
        )
        
        # Wait for both tasks to complete
        results = await asyncio.gather(
            alphamissense_task, clinvar_task, return_exceptions=True
        )
        
        # Handle exceptions
        alphamissense_result, clinvar_result = results
        
        if isinstance(alphamissense_result, Exception):
            alphamissense_result = AlphaMissenseResult(
                variant_id=variant_id,
                gene=gene,
                protein_change=protein_change,
                warnings=[f"AlphaMissense error: {str(alphamissense_result)}"]
            )
        
        if isinstance(clinvar_result, Exception):
            clinvar_result = ClinVarAnnotation(
                variant_id=variant_id,
                gene=gene,
                protein_change=protein_change,
                clinical_significance="Uncertain_Significance",
                warnings=[f"ClinVar error: {str(clinvar_result)}"]
            )
        else:
            # If clinical_significance is None or empty, set to "Uncertain_Significance"
            if getattr(clinvar_result, "clinical_significance", None) in (None, ""):
                object.__setattr__(clinvar_result, "clinical_significance", "Uncertain_Significance")
        # Ensure we have proper types
        assert isinstance(alphamissense_result, AlphaMissenseResult)
        assert isinstance(clinvar_result, ClinVarAnnotation)
        
        return VariantImpactResult(
            variant_id=variant_id,
            gene=gene,
            protein_change=protein_change,
            alphamissense=alphamissense_result,
            clinvar=clinvar_result
        )
    
    async def _analyze_alphamissense_async(self, variant_id: str, gene: str, 
                                          protein_change: str, vep_transcript_ids: List[str], 
                                          variant_data: Optional[Dict] = None) -> AlphaMissenseResult:
        """
        Analyze variant using AlphaMissense database asynchronously
        
        Args:
            variant_id: Variant identifier
            gene: Gene name
            protein_change: Protein change (e.g., "p.Arg175His")
            vep_transcript_ids: List of transcript IDs to search
            variant_data: Optional variant data containing reference and alternate alleles
            
        Returns:
            AlphaMissenseResult with analysis results
        """
        # Check cache first
        cached_result = self.cache.get_cached_result(variant_data or {})
        if cached_result:
            return cached_result
        
        # If in cache-only mode and no cached result, return default
        if self.cache_only_mode:
            logger.warning(f"No cached result found for {variant_id} in cache-only mode")
            return AlphaMissenseResult(
                variant_id=variant_id,
                gene=gene,
                protein_change=protein_change,
                average_pathogenicity_score=0.5,
                max_occurring_class=None,
                confidence="none",
                warnings=[f"No cached result available for {variant_id} in cache-only mode"]
            )
        
        # Run database query in thread pool to avoid blocking
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, self._analyze_alphamissense_sync, variant_id, gene, protein_change, vep_transcript_ids, variant_data
        )
        
        # Cache the result (but not if confidence is "none")
        if variant_data and result.confidence != "none":
            self.cache.cache_result(variant_data, result)
        
        return result
    
    def _calculate_confidence(self, filtering_method: str, num_matches: int) -> str:
        """Calculate confidence based on filtering method and match count"""
        if num_matches == 0:
            return "none"
        if filtering_method == "exact":
            return "high"
        elif filtering_method == "ref_alt" and num_matches <= 10:
            return "medium"
        elif filtering_method == "ref_alt" and num_matches > 10:
            return "low"
        elif filtering_method == "pathogenicity":
            return "low"
        else:
            return "low"

    def _analyze_alphamissense_sync(self, variant_id: str, gene: str, 
                                   protein_change: str, vep_transcript_ids: List[str], 
                                   variant_data: Optional[Dict] = None) -> AlphaMissenseResult:
        """
        Synchronous AlphaMissense analysis (runs in thread pool)
        
        Args:
            variant_id: Variant identifier
            gene: Gene name
            protein_change: Protein change (e.g., "p.Arg175His")
            vep_transcript_ids: List of transcript IDs to search
            variant_data: Optional variant data containing reference and alternate alleles
            
        Returns:
            AlphaMissenseResult with analysis results
        """
        result = AlphaMissenseResult(
            variant_id=variant_id,
            gene=gene,
            protein_change=protein_change
        )
        
        filtering_method = "none"
        matches = []
        try:
            # Extract position number from protein change
            position_match = re.search(r'(\d+)', protein_change)
            if not position_match:
                result.warnings.append(f"Could not extract position from protein change: {protein_change}")
                # If nothing works, default the pathogenic score to midway
                result.average_pathogenicity_score = 0.5
                result.confidence = "none"
                result.warnings.append("Defaulted pathogenicity score to 0.5 due to inability to extract position.")
                return result
            
            position_number = position_match.group(1)
            
            # Connect to database
            conn = sqlite3.connect(self.alphamissense_db_path)
            cursor = conn.cursor()
            
            # Search for variants in the specified transcripts with matching position
            query = """
                SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                FROM alphamissense 
                WHERE transcript_id IN ({})
                AND protein_variant LIKE ?
                ORDER BY am_pathogenicity DESC
            """.format(','.join(['?' for _ in vep_transcript_ids]))
            
            # Create pattern for position matching
            position_pattern = f"%{position_number}%"
            
            # Execute query
            params = vep_transcript_ids + [position_pattern]
            cursor.execute(query, params)
            matches = cursor.fetchall()
            if matches and len(matches) < 10:
                filtering_method = "exact"
            
            # If too many matches, try to filter by exact protein_variant
            if matches and len(matches) > 10:
                # Try to match the full protein_change string (e.g., 'p.Glu23fs')
                strict_query = """
                    SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                    FROM alphamissense 
                    WHERE transcript_id IN ({})
                    AND protein_variant = ?
                    ORDER BY am_pathogenicity DESC
                """.format(','.join(['?' for _ in vep_transcript_ids]))
                strict_params = vep_transcript_ids + [protein_change]
                cursor.execute(strict_query, strict_params)
                strict_matches = cursor.fetchall()
                if strict_matches:
                    matches = strict_matches
                    filtering_method = "exact"
                else:
                    # If strict matching fails, try filtering by reference and alternate alleles
                    if variant_data and 'reference' in variant_data and 'alternate' in variant_data:
                        ref_allele = variant_data['reference']
                        alt_allele = variant_data['alternate']
                        
                        ref_alt_query = """
                            SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                            FROM alphamissense 
                            WHERE transcript_id IN ({})
                            AND protein_variant LIKE ?
                            AND REF = ?
                            AND ALT = ?
                            ORDER BY am_pathogenicity DESC
                        """.format(','.join(['?' for _ in vep_transcript_ids]))
                        ref_alt_params = vep_transcript_ids + [position_pattern, ref_allele, alt_allele]
                        cursor.execute(ref_alt_query, ref_alt_params)
                        ref_alt_matches = cursor.fetchall()
                        if ref_alt_matches:
                            if len(ref_alt_matches) > 10:
                                result.warnings.append(f"REF/ALT filtering successful but too many matches ({len(ref_alt_matches)}), applying pathogenicity filtering")
                                # Filter by pathogenicity score on REF/ALT matches
                                pathogenicity_query = """
                                    SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                                    FROM alphamissense 
                                    WHERE transcript_id IN ({})
                                    AND protein_variant LIKE ?
                                    AND REF = ?
                                    AND ALT = ?
                                    ORDER BY am_pathogenicity DESC
                                    LIMIT 1
                                """.format(','.join(['?' for _ in vep_transcript_ids]))
                                pathogenicity_params = vep_transcript_ids + [position_pattern, ref_allele, alt_allele]
                                cursor.execute(pathogenicity_query, pathogenicity_params)
                                pathogenicity_matches = cursor.fetchall()
                                if pathogenicity_matches:
                                    matches = pathogenicity_matches
                                    filtering_method = "pathogenicity"
                                    result.warnings.append(f"Filtered to {len(matches)} high-pathogenicity matches (score > 0.8)")
                                else:
                                    matches = ref_alt_matches
                                    filtering_method = "ref_alt"
                                    result.warnings.append(f"Pathogenicity filtering failed, using REF/ALT matches: {len(matches)} matches")
                            else:
                                matches = ref_alt_matches
                                filtering_method = "ref_alt"
                                result.warnings.append(f"Filtered to {len(matches)} matches by reference/alternate alleles")
                        else:
                            # If reference/alternate filtering fails, try filtering by pathogenicity score
                            pathogenicity_query = """
                                SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                                FROM alphamissense 
                                WHERE transcript_id IN ({})
                                AND protein_variant LIKE ?
                                AND am_pathogenicity > 0.8
                                ORDER BY am_pathogenicity DESC
                                LIMIT 10
                            """.format(','.join(['?' for _ in vep_transcript_ids]))
                            pathogenicity_params = vep_transcript_ids + [position_pattern]
                            cursor.execute(pathogenicity_query, pathogenicity_params)
                            pathogenicity_matches = cursor.fetchall()
                            if pathogenicity_matches:
                                matches = pathogenicity_matches
                                filtering_method = "pathogenicity"
                                result.warnings.append(f"Filtered to {len(matches)} high-pathogenicity matches (score > 0.8)")
                            else:
                                result.warnings.append(f"Strict filtering by protein_variant and pathogenicity gave no results; using broader position-based matches.")
                    else:
                        # If no variant data available, try filtering by pathogenicity score
                        pathogenicity_query = """
                            SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                            FROM alphamissense 
                            WHERE transcript_id IN ({})
                            AND protein_variant LIKE ?
                            AND am_pathogenicity > 0.8
                            ORDER BY am_pathogenicity DESC
                            LIMIT 10
                        """.format(','.join(['?' for _ in vep_transcript_ids]))
                        pathogenicity_params = vep_transcript_ids + [position_pattern]
                        cursor.execute(pathogenicity_query, pathogenicity_params)
                        pathogenicity_matches = cursor.fetchall()
                        if pathogenicity_matches:
                            matches = pathogenicity_matches
                            filtering_method = "pathogenicity"
                            result.warnings.append(f"Filtered to {len(matches)} high-pathogenicity matches (score > 0.8)")
                        else:
                            result.warnings.append(f"Strict filtering by protein_variant and pathogenicity gave no results; using broader position-based matches.")
            
            if matches:
                # Extract pathogenicity scores and classes
                pathogenicity_scores = [match[0] for match in matches if match[0] is not None]
                classes = [match[1] for match in matches if match[1] is not None]
                matching_transcripts = list(set([match[3] for match in matches]))
                
                if pathogenicity_scores:
                    result.average_pathogenicity_score = statistics.mean(pathogenicity_scores)
                
                if classes:
                    # Find most common class
                    class_counter = Counter(classes)
                    result.max_occurring_class = class_counter.most_common(1)[0][0]
                
                result.matching_transcripts = matching_transcripts
                result.total_matches = len(matches)
                
    
            else:
                result.warnings.append(f"No AlphaMissense matches found for {variant_id}")
    
                # If nothing works, default the pathogenic score to midway
                result.average_pathogenicity_score = 0.5
                result.confidence = "none"
                result.warnings.append("Defaulted pathogenicity score to 0.5 due to no AlphaMissense matches found.")
            
            conn.close()
        
        except Exception as e:
            error_msg = f"Error analyzing AlphaMissense for {variant_id}: {str(e)}"
            result.warnings.append(error_msg)
            logger.error(error_msg)
            # If nothing works, default the pathogenic score to midway
            result.average_pathogenicity_score = 0.5
            result.confidence = "none"
            result.warnings.append("Defaulted pathogenicity score to 0.5 due to AlphaMissense error.")

        # Calculate confidence if not already set (avoid overwriting if set above)
        if not hasattr(result, "confidence") or result.confidence is None:
            result.confidence = self._calculate_confidence(filtering_method, len(matches))
        return result
    
    async def _analyze_clinvar_async(self, variant: Dict) -> ClinVarAnnotation:
        """
        Analyze variant using ClinVar asynchronously
        
        Args:
            variant: Variant dictionary from input
            
        Returns:
            ClinVarAnnotation with analysis results
        """
        # Run ClinVar analysis in thread pool to avoid blocking
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, self.clinvar_annotator.annotate_variant, variant
        )
        # Default clinical_significance to "Uncertain_Significance" if not present or None/empty
        if getattr(result, "clinical_significance", None) in (None, ""):
            result.clinical_significance = "Uncertain_Significance"
        return result


def main():
    """Main function for testing"""
    
    parser = argparse.ArgumentParser(description='Parallel Variant Impact Analyzer')
    parser.add_argument('--input', '-i', default="processed_input.json", 
                       help='Input file path (default: processed_input.json)')
    parser.add_argument('--output', '-o', default="variant_impact_results.json", 
                       help='Output file path (default: variant_impact_results.json)')
    parser.add_argument('--db-path', default="cache/alphamissense/alphamissense_hg38.db",
                       help='AlphaMissense database path (default: cache/alphamissense/alphamissense_hg38.db)')
    parser.add_argument('--workers', '-w', type=int, default=None,
                       help='Number of parallel workers (default: CPU count)')
    parser.add_argument('--cache-dir', default="cache",
                       help='Cache directory (default: cache)')
    
    args = parser.parse_args()
    
    analyzer = VariantImpactAnalyzer(
        alphamissense_db_path=args.db_path,
        max_workers=args.workers,
        cache_dir=args.cache_dir
    )
    
    # Analyze variants
    logger.info(f"Analyzing variants from: {args.input}")
    start_time = time.time()
    results = analyzer.analyze_variants(args.input)
    end_time = time.time()
    
    # Save results
    analyzer.save_results(results, args.output)
    
    # Print summary
    logger.info(f"Analysis complete in {end_time - start_time:.2f} seconds")
    
    # Print summary statistics
    variants = results.get('missense_variants', [])
    variants_with_scores = [v for v in variants if v.get('pathogenicity_score')]
    
    logger.info(f"Total variants processed: {len(variants)}")
    logger.info(f"Variants with pathogenicity scores: {len(variants_with_scores)}")
    
    if variants_with_scores:
        scores = [v['pathogenicity_score'] for v in variants_with_scores]
        logger.info(f"Average pathogenicity score: {sum(scores)/len(scores):.4f}")
    
    # Count genes
    genes = set(v['gene'] for v in variants)
    logger.info(f"Unique genes: {len(genes)}")
    
    # Show cache-only mode info if applicable
    if analyzer.cache_only_mode:
        logger.info("Analysis completed in cache-only mode (no AlphaMissense database available)")


if __name__ == "__main__":
    main() 