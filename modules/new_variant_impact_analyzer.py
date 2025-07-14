#!/usr/bin/env python3
"""
New Variant Impact Analyzer

This module processes variants in parallel using:
1. AlphaMissense database for pathogenicity scores
2. ClinVar database for clinical annotations

Features:
- Parallel processing of variants
- Fast AlphaMissense database queries
- ClinVar annotation lookup
- Comprehensive result aggregation
"""

import os
import json
import logging
import sqlite3
import re
import pandas as pd
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class AlphaMissenseResult:
    """Result from AlphaMissense database query"""
    variant_id: str
    transcript_ids: List[str]
    protein_change: str
    average_score: Optional[float] = None
    scores: List[float] = field(default_factory=list)
    matches: int = 0
    warnings: List[str] = field(default_factory=list)


@dataclass
class ClinVarResult:
    """Result from ClinVar database query"""
    variant_id: str
    gene: str
    protein_change: str
    annotations: List[Dict] = field(default_factory=list)
    significance: Optional[str] = None
    variation_id: Optional[str] = None
    warnings: List[str] = field(default_factory=list)


@dataclass
class VariantImpactResult:
    """Combined result for a variant"""
    variant_id: str
    gene: str
    protein_change: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    transcript_ids: List[str]
    
    # AlphaMissense results
    alphamissense_score: Optional[float] = None
    alphamissense_matches: int = 0
    
    # ClinVar results
    clinvar_annotations: List[Dict] = field(default_factory=list)
    clinvar_significance: Optional[str] = None
    
    # Combined data
    warnings: List[str] = field(default_factory=list)


class AlphaMissenseProcessor:
    """Process AlphaMissense database queries"""
    
    def __init__(self, db_path: str = "cache/alphamissense/alphamissense.db"):
        self.db_path = db_path
        self.conn = None
        self._connect_db()
    
    def _connect_db(self):
        """Connect to AlphaMissense database"""
        # Try multiple possible database paths
        possible_paths = [
            "cache/alphamissense/alphamissense.db",
            "cache/alphamissense/alphamissense_hg38.db",
            "alphamissense.db"
        ]
        
        for db_path in possible_paths:
            try:
                if os.path.exists(db_path):
                    self.db_path = db_path
                    self.conn = sqlite3.connect(db_path)
                    logger.info(f"Connected to AlphaMissense database: {db_path}")
                    return
            except Exception as e:
                logger.warning(f"Could not connect to {db_path}: {e}")
                continue
        
        # If no database found, create a minimal one for testing
        logger.warning("No AlphaMissense database found. Creating minimal test database.")
        self.db_path = "cache/alphamissense/alphamissense_test.db"
        self.conn = sqlite3.connect(self.db_path)
        self._create_test_db()
    
    def _create_test_db(self):
        """Create a minimal test database for development"""
        try:
            if self.conn is None:
                logger.error("No database connection available")
                return
                
            cursor = self.conn.cursor()
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS alphamissense (
                    gene TEXT,
                    transcript TEXT,
                    protein_variant TEXT,
                    am_pathogenicity REAL
                )
            """)
            
            # Insert some test data
            test_data = [
                ('TP53', 'ENST00000269305', 'R175H', 0.85),
                ('BRCA1', 'ENST00000352993', 'C61G', 0.72),
                ('CFTR', 'ENST00000003084', 'F508del', 0.91)
            ]
            
            cursor.executemany("""
                INSERT INTO alphamissense (gene, transcript, protein_variant, am_pathogenicity)
                VALUES (?, ?, ?, ?)
            """, test_data)
            
            self.conn.commit()
            logger.info("Created test AlphaMissense database with sample data")
            
        except Exception as e:
            logger.error(f"Error creating test database: {e}")
    
    def process_variant(self, variant: Dict) -> AlphaMissenseResult:
        """Process a single variant through AlphaMissense database"""
        variant_id = variant.get('id', 'unknown')
        gene = variant.get('gene', '')
        protein_change = variant.get('protein_change', '')
        transcript_ids = variant.get('vep_transcript_ids', [])
        
        result = AlphaMissenseResult(
            variant_id=variant_id,
            transcript_ids=transcript_ids,
            protein_change=protein_change
        )
        
        try:
            # Extract position number from protein change (e.g., "p.Arg175His" -> "175")
            position_number = self._extract_position_number(protein_change)
            if not position_number:
                result.warnings.append(f"Could not extract position number from protein change: {protein_change}")
                return result
            
            # Query AlphaMissense database
            scores = self._query_alphamissense(gene, position_number, transcript_ids)
            
            if scores:
                result.scores = scores
                result.average_score = sum(scores) / len(scores)
                result.matches = len(scores)
            else:
                result.warnings.append(f"No AlphaMissense matches found for {gene} position {position_number}")
                
        except Exception as e:
            result.warnings.append(f"Error processing AlphaMissense query: {e}")
            logger.error(f"AlphaMissense processing error for {variant_id}: {e}")
        
        return result
    
    def _extract_position_number(self, protein_change: str) -> Optional[str]:
        """Extract position number from protein change string"""
        # Pattern to match position number in protein change
        # Examples: "p.Arg175His" -> "175", "p.Cys61Gly" -> "61"
        pattern = r'p\.[A-Za-z]+(\d+)[A-Za-z]+'
        match = re.search(pattern, protein_change)
        if match:
            return match.group(1)
        return None
    
    def _query_alphamissense(self, gene: str, position_number: str, transcript_ids: List[str]) -> List[float]:
        """Query AlphaMissense database for scores"""
        if not self.conn:
            return []
        
        try:
            # Build query to search for variants with the position number
            # Search in protein_variant column for the position number
            query = """
            SELECT am_pathogenicity 
            FROM alphamissense 
            WHERE gene = ? 
            AND protein_variant LIKE ?
            """
            
            # Add transcript filter if available
            if transcript_ids:
                transcript_placeholders = ','.join(['?' for _ in transcript_ids])
                query += f" AND transcript IN ({transcript_placeholders})"
            
            # Prepare parameters
            params = [gene, f"%{position_number}%"]
            if transcript_ids:
                params.extend(transcript_ids)
            
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            
            results = cursor.fetchall()
            scores = [row[0] for row in results if row[0] is not None]
            
            return scores
            
        except Exception as e:
            logger.error(f"Error querying AlphaMissense database: {e}")
            return []
    
    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()


class ClinVarProcessor:
    """Process ClinVar database queries"""
    
    def __init__(self, cache_dir: str = "cache/clinvar"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.cache_db_path = self.cache_dir / "clinvar_cache.db"
        self._init_cache_db()
        
        # Use existing ClinVar local lookup as fallback
        from clinvar_local_lookup import ClinVarLocalLookup
        self.clinvar_lookup = ClinVarLocalLookup()
        
        # Initialize API access
        self.clinvar_api_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.api_key = os.getenv('NCBI_API_KEY', '')
        self.cache = {}
    
    def _init_cache_db(self):
        """Initialize ClinVar cache database"""
        try:
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            # Create cache table if it doesn't exist
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS clinvar_cache (
                    variant_key TEXT PRIMARY KEY,
                    gene TEXT,
                    protein_change TEXT,
                    significance TEXT,
                    variation_id TEXT,
                    annotations TEXT,
                    timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            conn.commit()
            conn.close()
            logger.info(f"Initialized ClinVar cache database: {self.cache_db_path}")
            
        except Exception as e:
            logger.error(f"Error initializing ClinVar cache: {e}")
    
    def process_variant(self, variant: Dict) -> ClinVarResult:
        """Process a single variant through ClinVar database"""
        variant_id = variant.get('id', 'unknown')
        gene = variant.get('gene', '')
        protein_change = variant.get('protein_change', '')
        
        result = ClinVarResult(
            variant_id=variant_id,
            gene=gene,
            protein_change=protein_change
        )
        
        try:
            # Try local lookup first
            if gene and protein_change:
                local_result = self.clinvar_lookup.lookup(gene, protein_change)
                if local_result:
                    result.significance = local_result.get('significance')
                    result.variation_id = local_result.get('variation_id')
                    result.annotations = [{
                        'significance': local_result.get('significance'),
                        'variation_id': local_result.get('variation_id'),
                        'raw_data': local_result.get('raw', [])
                    }]
                    return result
            
            # Fallback to API search
            cache_key = self._create_cache_key(variant)
            if cache_key in self.cache:
                cached_data = self.cache[cache_key]
                result.significance = cached_data.get('clinvar_significance')
                result.variation_id = cached_data.get('variation_id')
                result.annotations = [{
                    'significance': cached_data.get('clinvar_significance'),
                    'variation_id': cached_data.get('variation_id'),
                    'review_status': cached_data.get('review_status', ''),
                    'condition': cached_data.get('condition', '')
                }]
                return result
            
            # Search ClinVar API
            variant_ids = self._search_clinvar_ids(variant)
            if not variant_ids:
                result.warnings.append(f"No ClinVar annotations found for {gene} {protein_change}")
                return result
            
            # Fetch and parse all records, prioritize Pathogenic
            best_result = None
            for variant_id in variant_ids:
                clinvar_data = self._get_clinvar_data(variant_id)
                if not clinvar_data:
                    continue
                parsed_data = self._parse_clinvar_data(clinvar_data)
                if parsed_data['clinvar_significance'] == 'Pathogenic':
                    self.cache[cache_key] = parsed_data
                    result.significance = parsed_data['clinvar_significance']
                    result.variation_id = variant_id
                    result.annotations = [{
                        'significance': parsed_data['clinvar_significance'],
                        'variation_id': variant_id,
                        'review_status': parsed_data.get('review_status', ''),
                        'condition': parsed_data.get('condition', ''),
                        'last_evaluated': parsed_data.get('last_evaluated', '')
                    }]
                    return result
                if not best_result:
                    best_result = parsed_data
                    best_variant_id = variant_id
            
            if best_result:
                self.cache[cache_key] = best_result
                result.significance = best_result['clinvar_significance']
                result.variation_id = best_variant_id
                result.annotations = [{
                    'significance': best_result['clinvar_significance'],
                    'variation_id': best_variant_id,
                    'review_status': best_result.get('review_status', ''),
                    'condition': best_result.get('condition', ''),
                    'last_evaluated': best_result.get('last_evaluated', '')
                }]
            else:
                result.warnings.append(f"No ClinVar annotations found for {gene} {protein_change}")
                
        except Exception as e:
            result.warnings.append(f"Error processing ClinVar query: {e}")
            logger.error(f"ClinVar processing error for {variant_id}: {e}")
        
        return result
    
    def _get_cached_result(self, gene: str, protein_change: str) -> Optional[Dict]:
        """Get cached ClinVar result"""
        try:
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                SELECT significance, variation_id, annotations 
                FROM clinvar_cache 
                WHERE gene = ? AND protein_change = ?
            """, (gene, protein_change))
            
            result = cursor.fetchone()
            conn.close()
            
            if result:
                return {
                    'significance': result[0],
                    'variation_id': result[1],
                    'annotations': json.loads(result[2]) if result[2] else []
                }
            
        except Exception as e:
            logger.error(f"Error getting cached ClinVar result: {e}")
        
        return None
    
    def _cache_result(self, gene: str, protein_change: str, data: Dict):
        """Cache ClinVar result"""
        try:
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                INSERT OR REPLACE INTO clinvar_cache 
                (gene, protein_change, significance, variation_id, annotations)
                VALUES (?, ?, ?, ?, ?)
            """, (
                gene, 
                protein_change, 
                data.get('significance'),
                data.get('variation_id'),
                json.dumps(data.get('annotations', []))
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logger.error(f"Error caching ClinVar result: {e}")
    
    def _create_cache_key(self, variant: Dict) -> str:
        """Create cache key for variant"""
        gene = variant.get('gene', '')
        protein_change = variant.get('protein_change', '')
        chrom = variant.get('chromosome', '')
        pos = variant.get('position', '')
        ref = variant.get('reference_allele', '')
        alt = variant.get('alternate_allele', '')
        
        return f"{gene}_{protein_change}_{chrom}_{pos}_{ref}_{alt}"

    def _search_clinvar_ids(self, variant: Dict) -> Optional[list]:
        """Return all ClinVar IDs for a variant search"""
        try:
            gene = variant.get('gene', '')
            protein_change = variant.get('protein_change', '')
            pos = variant.get('position', '')
            ref = variant.get('reference_allele', '')
            alt = variant.get('alternate_allele', '')
            if not gene:
                return None
            queries = []
            if protein_change and pos and ref and alt:
                queries.append(f"{gene}[gene] AND {protein_change} AND {pos}[chrpos] AND {ref}[ref] AND {alt}[alt]")
            if protein_change and pos:
                queries.append(f"{gene}[gene] AND {protein_change} AND {pos}[chrpos]")
            if protein_change:
                queries.append(f"{gene}[gene] AND {protein_change}")
            queries.append(gene)
            for query in queries:
                logger.info(f"ClinVar search query: {query}")
                search_url = f"{self.clinvar_api_base}/esearch.fcgi"
                params = {
                    'db': 'clinvar',
                    'term': query,
                    'retmode': 'json',
                    'retmax': 10
                }
                if self.api_key:
                    params['api_key'] = self.api_key
                response = requests.get(search_url, params=params, timeout=30)
                if response.status_code != 200:
                    logger.warning(f"ClinVar search failed: {response.status_code}")
                    continue
                data = response.json()
                if 'esearchresult' not in data or 'idlist' not in data['esearchresult']:
                    continue
                id_list = data['esearchresult']['idlist']
                logger.info(f"ClinVar search returned IDs: {id_list}")
                if id_list:
                    return id_list
            return None
        except Exception as e:
            logger.warning(f"Error searching ClinVar: {e}")
            return None
    
    def _get_clinvar_data(self, variant_id: str) -> Optional[Dict]:
        """Get ClinVar data for variant using NCBI E-utilities, use rettype=vcv"""
        try:
            # Fetch ClinVar record
            fetch_url = f"{self.clinvar_api_base}/efetch.fcgi"
            params = {
                'db': 'clinvar',
                'id': variant_id,
                'rettype': 'vcv',
                'retmode': 'xml'
            }
            if self.api_key:
                params['api_key'] = self.api_key
            response = requests.get(fetch_url, params=params, timeout=30)
            if response.status_code != 200:
                logger.warning(f"ClinVar fetch failed: {response.status_code}")
                return None
            xml_content = response.text
            # Log XML for benchmark variants
            if any(x in xml_content for x in ["Arg175His", "Phe508del", "Glu23fs", "185delAG"]):
                with open(f"clinvar_debug_{variant_id}.xml", "w", encoding="utf-8") as f:
                    f.write(xml_content)
                logger.info(f"Saved ClinVar XML for benchmark variant {variant_id} to clinvar_debug_{variant_id}.xml")
            return self._parse_clinvar_xml(xml_content)
        except Exception as e:
            logger.warning(f"Error fetching ClinVar data: {e}")
            return None
    
    def _parse_clinvar_xml(self, xml_content: str) -> Optional[Dict]:
        """Parse ClinVar XML response using ElementTree and prioritize Pathogenic if present"""
        try:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(xml_content)
            # Find all ClinicalSignificance elements
            significances = []
            for cs in root.iter('ClinicalSignificance'):
                for child in cs:
                    if child.tag == 'Description':
                        significances.append(child.text)
            # Prioritize Pathogenic if present
            significance = 'Uncertain_significance'
            if significances:
                if 'Pathogenic' in significances:
                    significance = 'Pathogenic'
                elif 'Likely_pathogenic' in significances:
                    significance = 'Likely_pathogenic'
                elif 'Uncertain_significance' in significances:
                    significance = 'Uncertain_significance'
                elif 'Likely_benign' in significances:
                    significance = 'Likely_benign'
                elif 'Benign' in significances:
                    significance = 'Benign'
                else:
                    significance = significances[0]
            # Extract review status
            review_status = ''
            for rs in root.iter('ReviewStatus'):
                review_status = rs.text
                break
            # Extract condition
            condition = ''
            for trait in root.iter('Trait'):
                for name in trait.iter('Name'):
                    condition = name.text
                    break
                if condition:
                    break
            # Extract last evaluated date
            last_evaluated = ''
            for date in root.iter('DateLastEvaluated'):
                last_evaluated = date.text
                break
            return {
                'significance': significance,
                'review_status': review_status,
                'condition': condition,
                'last_evaluated': last_evaluated
            }
        except Exception as e:
            logger.warning(f"Error parsing ClinVar XML: {e}")
            return None
    
    def _parse_clinvar_data(self, data: Dict) -> Dict:
        """Parse ClinVar data"""
        significance = data.get('significance', 'Uncertain_significance')
        
        # Convert ClinVar significance to score
        significance_scores = {
            'Pathogenic': 0.9,
            'Likely_pathogenic': 0.8,
            'Uncertain_significance': 0.5,
            'Likely_benign': 0.2,
            'Benign': 0.1,
            'Conflicting_interpretations_of_pathogenicity': 0.5
        }
        
        score = significance_scores.get(significance, 0.5)
        
        return {
            'clinvar_significance': significance,
            'clinvar_score': score,
            'review_status': data.get('review_status', ''),
            'condition': data.get('condition', ''),
            'last_evaluated': data.get('last_evaluated', ''),
            'variation_id': data.get('variation_id')
        }

    def _significance_to_score(self, significance: str) -> float:
        s = significance.lower()
        if 'pathogenic' in s:
            return 0.9
        elif 'likely_pathogenic' in s:
            return 0.8
        elif 'uncertain' in s:
            return 0.5
        elif 'likely_benign' in s:
            return 0.2
        elif 'benign' in s:
            return 0.1
        else:
            return 0.5

    def get_complete_clinvar_records(self, variant_ids: List[str]) -> List[Dict]:
        """Get complete ClinVar records for a list of variation IDs"""
        complete_records = []
        for variant_id in variant_ids:
            record = self._get_clinvar_data(variant_id)
            if record:
                record['variation_id'] = variant_id
                complete_records.append(record)
        return complete_records

    def _query_clinvar(self, gene: str, protein_change: str) -> Optional[Dict]:
        """Query ClinVar database for annotations"""
        # This is a placeholder implementation
        # In a real implementation, you would query the actual ClinVar database
        # For now, we'll return a basic structure
        
        # Try to find exact match or similar matches
        # This would typically involve querying a local ClinVar database
        # or using the ClinVar API
        
        return {
            'significance': 'Unknown',
            'variation_id': None,
            'annotations': []
        }


class NewVariantImpactAnalyzer:
    """Main analyzer class for processing variants in parallel"""
    
    def __init__(self, max_workers: int = 4):
        self.max_workers = max_workers
        self.alphamissense_processor = AlphaMissenseProcessor()
        self.clinvar_processor = ClinVarProcessor()
    
    def analyze_variants(self, input_file: str = "processed_input.json") -> Dict[str, Any]:
        """Analyze variants from processed input file"""
        try:
            # Load processed input
            with open(input_file, 'r') as f:
                data = json.load(f)
            
            variants = data.get('missense_variants', [])
            logger.info(f"Processing {len(variants)} variants")
            
            # Process variants in parallel
            results = self._process_variants_parallel(variants)
            
            # Combine results with original data
            output = {
                'input_data': data,
                'variant_results': results,
                'summary': self._generate_summary(results)
            }
            
            return output
            
        except Exception as e:
            logger.error(f"Error analyzing variants: {e}")
            raise
    
    def _process_variants_parallel(self, variants: List[Dict]) -> Dict[str, VariantImpactResult]:
        """Process variants in parallel"""
        results = {}
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all variant processing tasks
            future_to_variant = {
                executor.submit(self._process_single_variant, variant): variant 
                for variant in variants
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_variant):
                variant = future_to_variant[future]
                try:
                    result = future.result()
                    results[variant['id']] = result
                    logger.info(f"Completed processing variant: {variant['id']}")
                except Exception as e:
                    logger.error(f"Error processing variant {variant['id']}: {e}")
                    # Create error result
                    results[variant['id']] = VariantImpactResult(
                        variant_id=variant['id'],
                        gene=variant.get('gene', ''),
                        protein_change=variant.get('protein_change', ''),
                        chromosome=variant.get('chromosome', ''),
                        position=variant.get('position', 0),
                        reference=variant.get('reference', ''),
                        alternate=variant.get('alternate', ''),
                        transcript_ids=variant.get('vep_transcript_ids', []),
                        warnings=[f"Processing error: {e}"]
                    )
        
        return results
    
    def _process_single_variant(self, variant: Dict) -> VariantImpactResult:
        """Process a single variant with both AlphaMissense and ClinVar"""
        variant_id = variant.get('id', 'unknown')
        
        # Process AlphaMissense and ClinVar in parallel
        with ThreadPoolExecutor(max_workers=2) as executor:
            am_future = executor.submit(self.alphamissense_processor.process_variant, variant)
            cv_future = executor.submit(self.clinvar_processor.process_variant, variant)
            
            # Wait for both results
            am_result = am_future.result()
            cv_result = cv_future.result()
        
        # Combine results
        result = VariantImpactResult(
            variant_id=variant_id,
            gene=variant.get('gene', ''),
            protein_change=variant.get('protein_change', ''),
            chromosome=variant.get('chromosome', ''),
            position=variant.get('position', 0),
            reference=variant.get('reference', ''),
            alternate=variant.get('alternate', ''),
            transcript_ids=variant.get('vep_transcript_ids', []),
            alphamissense_score=am_result.average_score,
            alphamissense_matches=am_result.matches,
            clinvar_annotations=cv_result.annotations,
            clinvar_significance=cv_result.significance,
            warnings=am_result.warnings + cv_result.warnings
        )
        
        return result
    
    def _generate_summary(self, results: Dict[str, VariantImpactResult]) -> Dict[str, Any]:
        """Generate summary statistics"""
        total_variants = len(results)
        variants_with_am_scores = sum(1 for r in results.values() if r.alphamissense_score is not None)
        variants_with_cv_annotations = sum(1 for r in results.values() if r.clinvar_annotations)
        
        return {
            'total_variants': total_variants,
            'variants_with_alphamissense_scores': variants_with_am_scores,
            'variants_with_clinvar_annotations': variants_with_cv_annotations,
            'alphamissense_coverage': variants_with_am_scores / total_variants if total_variants > 0 else 0,
            'clinvar_coverage': variants_with_cv_annotations / total_variants if total_variants > 0 else 0
        }
    
    def save_results(self, results: Dict[str, Any], output_file: str = "new_variant_impact_results.json"):
        """Save results to JSON file"""
        try:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            logger.info(f"Results saved to: {output_file}")
        except Exception as e:
            logger.error(f"Error saving results: {e}")
            raise
    
    def close(self):
        """Clean up resources"""
        if hasattr(self, 'alphamissense_processor'):
            self.alphamissense_processor.close()


def main():
    """Main function for testing"""
    analyzer = NewVariantImpactAnalyzer()
    
    try:
        # Analyze variants
        results = analyzer.analyze_variants()
        
        # Save results
        analyzer.save_results(results)
        
        # Print summary
        summary = results['summary']
        print(f"\nAnalysis Summary:")
        print(f"Total variants: {summary['total_variants']}")
        print(f"AlphaMissense coverage: {summary['alphamissense_coverage']:.2%}")
        print(f"ClinVar coverage: {summary['clinvar_coverage']:.2%}")
        
    except Exception as e:
        logger.error(f"Error in main: {e}")
    finally:
        analyzer.close()


if __name__ == "__main__":
    # TEMPORARY: Test only the ClinVar lookup for all variants in processed_input.json
    import sys
    print("Testing ClinVar lookup only...")
    try:
        with open("../processed_input.json", "r") as f:
            data = json.load(f)
        variants = data.get("missense_variants", [])
        clinvar_proc = ClinVarProcessor()
        
        print(f"\n{'='*80}")
        print("CLINVAR ANNOTATION RESULTS")
        print(f"{'='*80}")
        
        for variant in variants:
            result = clinvar_proc.process_variant(variant)
            
            print(f"\n{'='*60}")
            print(f"VARIANT: {variant['id']}")
            print(f"{'='*60}")
            print(f"Gene: {variant['gene']}")
            print(f"Protein Change: {variant['protein_change']}")
            print(f"Chromosome: {variant.get('chromosome', 'N/A')}")
            print(f"Position: {variant.get('position', 'N/A')}")
            print(f"Reference: {variant.get('reference', 'N/A')}")
            print(f"Alternate: {variant.get('alternate', 'N/A')}")
            
            print(f"\nCLINVAR RESULTS:")
            print(f"  Significance: {result.significance}")
            print(f"  Variation ID: {result.variation_id}")
            
            # Get all ClinVar IDs found during search
            variant_ids = clinvar_proc._search_clinvar_ids(variant)
            if variant_ids:
                print(f"  All ClinVar IDs found: {variant_ids}")
                print(f"  Complete ClinVar Records:")
                complete_records = clinvar_proc.get_complete_clinvar_records(variant_ids)
                for i, record in enumerate(complete_records, 1):
                    print(f"    Record {i} (ID: {record.get('variation_id', 'N/A')}):")
                    print(f"      Significance: {record.get('significance', 'N/A')}")
                    print(f"      Review Status: {record.get('review_status', 'N/A')}")
                    print(f"      Condition: {record.get('condition', 'N/A')}")
                    print(f"      Last Evaluated: {record.get('last_evaluated', 'N/A')}")
            else:
                print(f"  All ClinVar IDs found: None")
            
            if result.annotations:
                print(f"  Selected Annotation:")
                for i, annotation in enumerate(result.annotations, 1):
                    print(f"    Annotation {i}:")
                    print(f"      Significance: {annotation.get('significance', 'N/A')}")
                    print(f"      Variation ID: {annotation.get('variation_id', 'N/A')}")
                    print(f"      Review Status: {annotation.get('review_status', 'N/A')}")
                    print(f"      Condition: {annotation.get('condition', 'N/A')}")
                    print(f"      Last Evaluated: {annotation.get('last_evaluated', 'N/A')}")
                    if 'raw_data' in annotation:
                        print(f"      Raw Data: {annotation['raw_data']}")
            else:
                print(f"  Annotations: None found")
            
            if result.warnings:
                print(f"  Warnings: {result.warnings}")
            else:
                print(f"  Warnings: None")
                
        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")
        print(f"Total variants processed: {len(variants)}")
        print(f"Variants with ClinVar annotations: {sum(1 for v in variants if ClinVarProcessor().process_variant(v).annotations)}")
        print(f"Variants without ClinVar annotations: {sum(1 for v in variants if not ClinVarProcessor().process_variant(v).annotations)}")
        
    except Exception as e:
        print(f"Error during ClinVar test: {e}", file=sys.stderr) 