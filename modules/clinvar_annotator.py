#!/usr/bin/env python3
"""
ClinVar Annotator

A comprehensive module for annotating variants using ClinVar database.
Supports multiple lookup strategies:
1. Local ClinVar database cache
2. NCBI ClinVar API
3. Local variant summary file
4. Fallback mechanisms

Example usage:
    annotator = ClinVarAnnotator()
    result = annotator.annotate_variant({
        "id": "var_001",
        "position": 7675088,
        "reference": "G",
        "alternate": "A",
        "gene": "TP53",
        "protein_change": "p.Arg175His",
        "chromosome": "17",
        "transcript_id": "ENST00000269305"
    })
"""

import os
import json
import logging
import requests
import sqlite3
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import time
import xml.etree.ElementTree as ET

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class ClinVarAnnotation:
    """ClinVar variant annotation result"""
    variant_id: str
    gene: str
    protein_change: str
    clinical_significance: Optional[str] = None
    review_status: Optional[str] = None
    condition: Optional[str] = None
    variation_id: Optional[str] = None
    last_evaluated: Optional[str] = None
    submitter: Optional[str] = None
    star_rating: Optional[int] = None
    accession: Optional[str] = None
    classification_type: Optional[str] = None  # "germline", "somatic", "oncogenicity", "main"
    warnings: List[str] = field(default_factory=list)
    source: str = "unknown"  # "local_cache", "api", "local_file", "fallback"


class ClinVarAnnotator:
    """Comprehensive ClinVar variant annotator"""
    
    def __init__(self, cache_dir: str = "cache/clinvar", 
                 summary_file: str = "cache/clinvar_variant_summary.txt.gz"):
        """
        Initialize ClinVar annotator
        
        Args:
            cache_dir: Directory for caching ClinVar data
            summary_file: Path to local ClinVar variant summary file
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.summary_file = Path(summary_file)
        
        # API configuration
        self.api_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.api_key = os.getenv('NCBI_API_KEY', '34bce239811b7a833431a406c588de697008')
        
        # Initialize components
        self._init_cache_db()
        self._init_local_lookup()
        
        # Known pathogenic variants cache
        self.known_variants = self._init_known_variants()
        
        logger.info("ClinVar Annotator initialized")
    
    def _init_cache_db(self):
        """Initialize SQLite cache database"""
        try:
            self.cache_db_path = self.cache_dir / "clinvar_cache.db"
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            # Check if table exists and has the new column
            cursor.execute("PRAGMA table_info(variant_annotations)")
            columns = [col[1] for col in cursor.fetchall()]
            
            if 'classification_type' not in columns:
                # Drop and recreate table with new schema
                cursor.execute("DROP TABLE IF EXISTS variant_annotations")
                logger.info("Recreating ClinVar cache database with updated schema")
            
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS variant_annotations (
                    variant_key TEXT PRIMARY KEY,
                    gene TEXT,
                    protein_change TEXT,
                    clinical_significance TEXT,
                    review_status TEXT,
                    condition TEXT,
                    variation_id TEXT,
                    last_evaluated TEXT,
                    submitter TEXT,
                    star_rating INTEGER,
                    accession TEXT,
                    classification_type TEXT,
                    source TEXT,
                    timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            conn.commit()
            conn.close()
            logger.info(f"Initialized ClinVar cache database: {self.cache_db_path}")
            
        except Exception as e:
            logger.error(f"Error initializing cache database: {e}")
    
    def _init_local_lookup(self):
        """Initialize local ClinVar lookup from summary file"""
        self.local_index = {}
        if self.summary_file.exists():
            try:
                # self._build_local_index()
                logger.info(f"Built local ClinVar index with {len(self.local_index)} entries")
            except Exception as e:
                logger.warning(f"Error building local index: {e}")
    
    def _build_local_index(self):
        """Build index from local ClinVar summary file"""
        # Check if index already exists and summary file hasn't changed
        index_cache_file = self.cache_dir / "local_index_cache.json"
        summary_mtime = self.summary_file.stat().st_mtime if self.summary_file.exists() else 0
        
        if index_cache_file.exists():
            try:
                with open(index_cache_file, 'r') as f:
                    cache_data = json.load(f)
                    if cache_data.get('summary_mtime', 0) == summary_mtime:
                        self.local_index = cache_data.get('index', {})
                        logger.info(f"Loaded cached local index with {len(self.local_index)} entries")
                        return
            except Exception as e:
                logger.warning(f"Error loading cached index: {e}")
        
        # Build index from scratch
        self.local_index = {}
        if self.summary_file.exists():
            with gzip.open(self.summary_file, 'rt', encoding='utf-8') as f:
                header = next(f).strip().split('\t')
                col_idx = {col: i for i, col in enumerate(header)}
                
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) < len(header):
                        continue
                    
                    gene = fields[col_idx.get('GeneSymbol', 0)]
                    protein = fields[col_idx.get('Name', 1)]
                    significance = fields[col_idx.get('ClinicalSignificance', 2)]
                    variation_id = fields[col_idx.get('VariationID', 3)]
                    
                    if gene and protein:
                        key = f"{gene}|{protein}"
                        self.local_index[key] = {
                            'significance': significance,
                            'variation_id': variation_id,
                            'raw': fields
                        }
            
            # Cache the built index
            try:
                cache_data = {
                    'summary_mtime': summary_mtime,
                    'index': self.local_index
                }
                with open(index_cache_file, 'w') as f:
                    json.dump(cache_data, f)
                logger.info(f"Built and cached local ClinVar index with {len(self.local_index)} entries")
            except Exception as e:
                logger.warning(f"Error caching local index: {e}")
        else:
            logger.warning(f"ClinVar summary file not found: {self.summary_file}")
    
    def _init_known_variants(self) -> Dict[str, Dict]:
        """Initialize known pathogenic variants cache"""
        return {
           
        }
    
    def annotate_variant(self, variant: Dict) -> ClinVarAnnotation:
        """
        Annotate a single variant using ClinVar
        
        Args:
            variant: Dictionary containing variant information
            
        Returns:
            ClinVarAnnotation object with annotation results
        """
        variant_id = variant.get('id', 'unknown')
        gene = variant.get('gene', '')
        protein_change = variant.get('protein_change', '')
        
        annotation = ClinVarAnnotation(
            variant_id=variant_id,
            gene=gene,
            protein_change=protein_change
        )
        
        try:
            # Strategy 1: Check known variants cache
            if gene and protein_change:
                known_key = f"{gene}|{protein_change}"
                if known_key in self.known_variants:
                    known_data = self.known_variants[known_key]
                    annotation.clinical_significance = known_data['clinical_significance']
                    annotation.review_status = known_data['review_status']
                    annotation.condition = known_data['condition']
                    annotation.variation_id = known_data['variation_id']
                    annotation.star_rating = known_data['star_rating']
                    annotation.accession = known_data['accession']
                    annotation.source = "known_variants"
                    return annotation
            
            # Strategy 2: Check local cache database
            cached_result = self._get_cached_annotation(variant)
            if cached_result:
                annotation.clinical_significance = cached_result['clinical_significance']
                annotation.review_status = cached_result['review_status']
                annotation.condition = cached_result['condition']
                annotation.variation_id = cached_result['variation_id']
                annotation.last_evaluated = cached_result['last_evaluated']
                annotation.submitter = cached_result['submitter']
                annotation.star_rating = cached_result['star_rating']
                annotation.accession = cached_result['accession']
                annotation.classification_type = cached_result.get('classification_type', 'main')
                annotation.source = cached_result['source']
                return annotation
            
            # Strategy 3: Check local summary file
            local_result = self._lookup_local(variant)
            if local_result:
                annotation.clinical_significance = local_result['significance']
                annotation.variation_id = local_result['variation_id']
                annotation.source = "local_file"
                self._cache_annotation(variant, annotation)
                return annotation
            
            # Strategy 4: Query ClinVar API
            try:
                api_result = self._query_clinvar_api(variant)
                if api_result:
                    annotation.clinical_significance = api_result['clinical_significance']
                    annotation.review_status = api_result['review_status']
                    annotation.condition = api_result['condition']
                    annotation.variation_id = api_result['variation_id']
                    annotation.last_evaluated = api_result['last_evaluated']
                    annotation.submitter = api_result['submitter']
                    annotation.star_rating = api_result['star_rating']
                    annotation.accession = api_result['accession']
                    annotation.classification_type = api_result.get('classification_type', 'main')
                    annotation.source = "api"
                    self._cache_annotation(variant, annotation)
                    return annotation
            except Exception as e:
                logger.warning(f"Error in API query for {variant_id}: {e}")
                annotation.warnings.append(f"API query failed: {e}")
            
            # Strategy 5: Fallback - no annotation found
            annotation.warnings.append(f"No ClinVar annotation found for {gene} {protein_change}")
            annotation.source = "no_annotation"
            
        except Exception as e:
            annotation.warnings.append(f"Error during ClinVar annotation: {e}")
            logger.error(f"ClinVar annotation error for {variant_id}: {e}")
        
        return annotation
    
    def _get_cached_annotation(self, variant: Dict) -> Optional[Dict]:
        """Get cached annotation from database"""
        try:
            cache_key = self._create_cache_key(variant)
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                SELECT clinical_significance, review_status, condition, variation_id,
                       last_evaluated, submitter, star_rating, accession, classification_type, source
                FROM variant_annotations 
                WHERE variant_key = ?
            """, (cache_key,))
            
            result = cursor.fetchone()
            conn.close()
            
            if result:
                return {
                    'clinical_significance': result[0],
                    'review_status': result[1],
                    'condition': result[2],
                    'variation_id': result[3],
                    'last_evaluated': result[4],
                    'submitter': result[5],
                    'star_rating': result[6],
                    'accession': result[7],
                    'classification_type': result[8],
                    'source': result[9]
                }
            
        except Exception as e:
            logger.error(f"Error getting cached annotation: {e}")
        
        return None
    
    def _lookup_local(self, variant: Dict) -> Optional[Dict]:
        """Lookup variant in local summary file"""
        gene = variant.get('gene', '')
        protein_change = variant.get('protein_change', '')
        
        if not gene or not protein_change:
            return None
        
        # Try exact match
        key = f"{gene}|{protein_change}"
        if key in self.local_index:
            return self.local_index[key]
        
        # Try partial match
        for k, v in self.local_index.items():
            if k.startswith(f"{gene}|") and protein_change in k:
                return v
        
        return None
    
    def _query_clinvar_api(self, variant: Dict) -> Optional[Dict]:
        """Query ClinVar API for variant information"""
        try:
            # Search for variant IDs
            variant_ids = self._search_clinvar_ids(variant)
            logger.info(f"Found variant IDs: {variant_ids}")
            if not variant_ids:
                return None
            
            # Get detailed information for each variant ID
            for variant_id in variant_ids:
                try:
                    result = self._get_clinvar_details(variant_id)
                    logger.info(f"Processing variant ID: {variant_id}")
                    if result and result.get('clinical_significance') != 'Cannot_annotate':
                        logger.info(f"Found valid result for variant ID {variant_id}: {result.get('clinical_significance')}")
                        return result
                    elif result:
                        logger.info(f"Found result but clinical significance is 'Cannot_annotate' for variant ID {variant_id}")
                except Exception as e:
                    logger.warning(f"Error processing variant ID {variant_id}: {e}")
                    continue  # Try next variant ID
            
        except Exception as e:
            logger.warning(f"Error querying ClinVar API: {e}")
        
        return None
    
    def _search_clinvar_ids(self, variant: Dict) -> List[str]:
        """Search ClinVar for variant IDs"""
        try:
            gene = variant.get('gene', '')
            protein_change = variant.get('protein_change', '')
            position = variant.get('position', '')
            reference = variant.get('reference', '')
            alternate = variant.get('alternate', '')
            
            if not gene:
                return []
            
            # Build search queries in order of specificity
            queries = []
            
            if protein_change and position:
                queries.append(f"{gene}[gene] AND {protein_change} AND {position}[chrpos]")
            
            if protein_change:
                queries.append(f"{gene}[gene] AND {protein_change}")
            
            queries.append(gene)
            
            for query in queries:
                logger.info(f"ClinVar search query: {query}")
                
                search_url = f"{self.api_base}/esearch.fcgi"
                params = {
                    'db': 'clinvar',
                    'term': query,
                    'retmode': 'json',
                    'retmax': 10
                }
                
                if self.api_key:
                    params['api_key'] = self.api_key
                
                response = requests.get(search_url, params=params, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    # Remove debug prints to reduce noise
                    # print(data)
                    # print("====")
                    if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                        id_list = data['esearchresult']['idlist']
                        if id_list:
                            logger.info(f"ClinVar search returned IDs: {id_list}")
                            return id_list
                    else:
                        logger.warning(f"Unexpected ClinVar search response format: {data}")
                else:
                    logger.warning(f"ClinVar search failed with status code: {response.status_code}")
                
                time.sleep(0.1)  # Rate limiting
            
        except Exception as e:
            logger.warning(f"Error searching ClinVar: {e}")
        
        return []
    
    def _get_clinvar_details(self, variant_id: str) -> Optional[Dict]:
        """Get detailed ClinVar information for a variant ID"""
        try:
            # Try esummary first (more reliable for getting basic info)
            summary_url = f"{self.api_base}/esummary.fcgi"
            summary_params = {
                'db': 'clinvar',
                'id': variant_id,
                'retmode': 'json'
            }
            
            if self.api_key:
                summary_params['api_key'] = self.api_key
            
            logger.info(f"Fetching ClinVar summary for variant ID: {variant_id}")
            summary_response = requests.get(summary_url, params=summary_params, timeout=30)
            
            if summary_response.status_code == 200:
                summary_data = summary_response.json()
                if 'result' in summary_data and variant_id in summary_data['result']:
                    variant_data = summary_data['result'][variant_id]
                    logger.info(f"Successfully fetched ClinVar summary for variant ID: {variant_id}")
                    return self._parse_clinvar_summary(variant_data)
            
            # Fallback to efetch if esummary fails
            fetch_url = f"{self.api_base}/efetch.fcgi"
            fetch_params = {
                'db': 'clinvar',
                'id': variant_id,
                'rettype': 'variation',
                'retmode': 'xml'
            }
            
            if self.api_key:
                fetch_params['api_key'] = self.api_key
            
            logger.info(f"Trying efetch for variant ID: {variant_id}")
            fetch_response = requests.get(fetch_url, params=fetch_params, timeout=30)
            
            if fetch_response.status_code == 200:
                logger.info(f"Successfully fetched ClinVar XML for variant ID: {variant_id}")
                return self._parse_clinvar_xml(fetch_response.text)
            else:
                logger.warning(f"ClinVar fetch failed for variant ID {variant_id} with status code: {fetch_response.status_code}")
                logger.warning(f"Response text: {fetch_response.text[:500]}")
            
        except Exception as e:
            logger.warning(f"Error fetching ClinVar details for variant ID {variant_id}: {e}")
        
        return None
    
    def _parse_clinvar_xml(self, xml_content: str) -> Optional[Dict]:
        """Parse ClinVar XML response with support for new germline classification model"""
        try:
            # Check if XML is empty or malformed
            if not xml_content or xml_content.strip() == '':
                logger.warning("Empty XML response from ClinVar API")
                return None
            
            # Check for empty ClinVar result set
            if '<ClinVarResult-Set><set/></ClinVarResult-Set>' in xml_content:
                logger.warning("ClinVar returned empty result set - no data found for this variant ID")
                return None
            
            root = ET.fromstring(xml_content)
            
            # Check if root is empty
            if len(root) == 0:
                logger.warning("Empty ClinVar XML response")
                return None
            
            # Check if we have any ClinVar records
            clinvar_records = root.findall('.//ClinVarSet')
            if not clinvar_records:
                logger.warning("No ClinVarSet records found in XML response")
                return None
            
            # Initialize variables
            clinical_significance = 'Cannot_annotate'
            review_status = ''
            condition = ''
            last_evaluated = ''
            submitter = ''
            star_rating = 1
            accession = ''
            
            # Extract germline classification (new model since Jan 2024)
            germline_found = False
            classification_type = 'main'  # Default fallback
            for cs in root.iter('ClinicalSignificance'):
                # First try to find germline-specific classification
                germline_elem = cs.find('GermlineClassification')
                if germline_elem is not None:
                    desc_elem = germline_elem.find('Description')
                    if desc_elem is not None and desc_elem.text:
                        clinical_significance = desc_elem.text
                        germline_found = True
                        classification_type = 'germline'
                        logger.info(f"Found germline classification: {clinical_significance}")
                        
                        # Get germline review status
                        germline_review = germline_elem.find('ReviewStatus')
                        if germline_review is not None and germline_review.text:
                            review_status = germline_review.text
                        break
                
                # If no germline classification found, try somatic clinical impact
                if not germline_found:
                    somatic_elem = cs.find('SomaticClinicalImpact')
                    if somatic_elem is not None:
                        desc_elem = somatic_elem.find('Description')
                        if desc_elem is not None and desc_elem.text:
                            clinical_significance = desc_elem.text
                            classification_type = 'somatic'
                            logger.info(f"Found somatic clinical impact: {clinical_significance}")
                            
                            # Get somatic review status
                            somatic_review = somatic_elem.find('ReviewStatus')
                            if somatic_review is not None and somatic_review.text:
                                review_status = somatic_review.text
                            break
            
            # Fallback: extract from main ClinicalSignificance (old model)
            if not germline_found and clinical_significance == 'Cannot_annotate':
                for cs in root.iter('ClinicalSignificance'):
                    desc_elem = cs.find('Description')
                    if desc_elem is not None and desc_elem.text:
                        clinical_significance = desc_elem.text
                        logger.info(f"Found main clinical significance: {clinical_significance}")
                        
                        # Get main review status
                        main_review = cs.find('ReviewStatus')
                        if main_review is not None and main_review.text:
                            review_status = main_review.text
                        break
            
            # Prioritize Pathogenic if multiple significances found
            if clinical_significance != 'Cannot_annotate':
                if 'Pathogenic' in clinical_significance:
                    clinical_significance = 'Pathogenic'
                elif 'Likely_pathogenic' in clinical_significance:
                    clinical_significance = 'Likely_pathogenic'
                elif 'Uncertain_significance' in clinical_significance:
                    clinical_significance = 'Uncertain_significance'
                elif 'Likely_benign' in clinical_significance:
                    clinical_significance = 'Likely_benign'
                elif 'Benign' in clinical_significance:
                    clinical_significance = 'Benign'
            
            # Extract condition
            for trait in root.iter('Trait'):
                for name in trait.iter('Name'):
                    condition = name.text
                    break
                if condition:
                    break
            
            # Extract last evaluated date
            for date in root.iter('DateLastEvaluated'):
                last_evaluated = date.text
                break
            
            # Extract submitter
            for sub in root.iter('Submitter'):
                for name in sub.iter('Name'):
                    submitter = name.text
                    break
                if submitter:
                    break
            
            # Extract accession
            for acc in root.iter('Accession'):
                accession = acc.text
                break
            
            # Determine star rating based on review status
            if review_status:
                if 'practice guideline' in review_status.lower():
                    star_rating = 4
                elif 'expert panel' in review_status.lower():
                    star_rating = 3
                elif 'multiple submitters' in review_status.lower():
                    star_rating = 2
                elif 'single submitter' in review_status.lower():
                    star_rating = 1
                else:
                    star_rating = 0
            
            return {
                'clinical_significance': clinical_significance,
                'review_status': review_status,
                'condition': condition,
                'last_evaluated': last_evaluated,
                'submitter': submitter,
                'star_rating': star_rating,
                'accession': accession,
                'classification_type': classification_type
            }
            
        except Exception as e:
            logger.warning(f"Error parsing ClinVar XML: {e}")
            return None
    
    def _parse_clinvar_summary(self, variant_data: Dict) -> Optional[Dict]:
        """Parse ClinVar summary data from JSON response"""
        try:
            logger.info(f"ClinVar summary data structure: {list(variant_data.keys())}")

            # Use germline_classification if present
            clinical_significance = variant_data.get('germline_classification', None)
            classification_type = 'germline'
            if not clinical_significance or clinical_significance == '-' or clinical_significance is None:
                # Fallback to somatic
                clinical_significance = variant_data.get('clinical_impact_classification', None)
                classification_type = 'somatic'
            if not clinical_significance or clinical_significance == '-' or clinical_significance is None:
                # Fallback to oncogenicity
                clinical_significance = variant_data.get('oncogenicity_classification', None)
                classification_type = 'oncogenicity'
            if not clinical_significance or clinical_significance == '-' or clinical_significance is None:
                clinical_significance = 'Cannot_annotate'
                classification_type = 'main'

            review_status = variant_data.get('review_status', '')
            condition = variant_data.get('condition', '')
            last_evaluated = variant_data.get('last_evaluated', '')
            submitter = variant_data.get('submitter', '')
            accession = variant_data.get('accession', '')

            # Determine star rating based on review status
            star_rating = 1  # Default
            if review_status:
                if 'practice guideline' in review_status.lower():
                    star_rating = 4
                elif 'expert panel' in review_status.lower():
                    star_rating = 3
                elif 'multiple submitters' in review_status.lower():
                    star_rating = 2
                elif 'single submitter' in review_status.lower():
                    star_rating = 1
                else:
                    star_rating = 0

            # Prioritize Pathogenic if present
            if clinical_significance != 'Cannot_annotate':
                if 'Pathogenic' in clinical_significance:
                    clinical_significance = 'Pathogenic'
                elif 'Likely_pathogenic' in clinical_significance:
                    clinical_significance = 'Likely_pathogenic'
                elif 'Uncertain_significance' in clinical_significance:
                    clinical_significance = 'Uncertain_significance'
                elif 'Likely_benign' in clinical_significance:
                    clinical_significance = 'Likely_benign'
                elif 'Benign' in clinical_significance:
                    clinical_significance = 'Benign'

            return {
                'clinical_significance': clinical_significance,
                'review_status': review_status,
                'condition': condition,
                'last_evaluated': last_evaluated,
                'submitter': submitter,
                'star_rating': star_rating,
                'accession': accession,
                'classification_type': classification_type
            }
        except Exception as e:
            logger.warning(f"Error parsing ClinVar summary: {e}")
            return None
    
    def _create_cache_key(self, variant: Dict) -> str:
        """Create cache key for variant"""
        gene = variant.get('gene', '')
        protein_change = variant.get('protein_change', '')
        chromosome = variant.get('chromosome', '')
        position = variant.get('position', '')
        reference = variant.get('reference', '')
        alternate = variant.get('alternate', '')
        
        return f"{gene}_{protein_change}_{chromosome}_{position}_{reference}_{alternate}"
    
    def _cache_annotation(self, variant: Dict, annotation: ClinVarAnnotation):
        """Cache annotation in database"""
        try:
            cache_key = self._create_cache_key(variant)
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                INSERT OR REPLACE INTO variant_annotations 
                (variant_key, gene, protein_change, clinical_significance, review_status,
                 condition, variation_id, last_evaluated, submitter, star_rating, accession, classification_type, source)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                cache_key,
                annotation.gene,
                annotation.protein_change,
                annotation.clinical_significance,
                annotation.review_status,
                annotation.condition,
                annotation.variation_id,
                annotation.last_evaluated,
                annotation.submitter,
                annotation.star_rating,
                annotation.accession,
                annotation.classification_type,
                annotation.source
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logger.error(f"Error caching annotation: {e}")
    
    def annotate_variants(self, variants: List[Dict]) -> Dict[str, ClinVarAnnotation]:
        """
        Annotate multiple variants
        
        Args:
            variants: List of variant dictionaries
            
        Returns:
            Dictionary mapping variant IDs to ClinVarAnnotation objects
        """
        results = {}
        
        for variant in variants:
            variant_id = variant.get('id', 'unknown')
            results[variant_id] = self.annotate_variant(variant)
        
        return results


def main():
    """Test the ClinVar annotator with the provided example"""
    
    # Example variant from the user
    test_variant = {
        "id": "var_001",
        "position": 7675088,
        "reference": "G",
        "alternate": "A",
        "gene": "TP53",
        "protein_change": "p.Arg175His",
        "chromosome": "17",
        "transcript": "",
        "transcript_id": "ENST00000269305",
        "vep_transcript_ids": [
            "ENST00000269305",
            "ENST00000359597",
            "ENST00000413465",
            "ENST00000420246",
            "ENST00000445888",
            "ENST00000455263",
            "ENST00000503591",
            "ENST00000504290",
            "ENST00000504937",
            "ENST00000505014",
            "ENST00000508793",
            "ENST00000509690",
            "ENST00000510385",
            "ENST00000514944",
            "ENST00000574684",
            "ENST00000576024",
            "ENST00000604348",
            "ENST00000610292",
            "ENST00000610538",
            "ENST00000610623",
            "ENST00000618944",
            "ENST00000619186",
            "ENST00000619485",
            "ENST00000620739",
            "ENST00000622645",
            "ENST00000635293",
            "ENST00000714356",
            "ENST00000714357",
            "ENST00000714358",
            "ENST00000714359",
            "ENST00000714408",
            "ENST00000714409"
        ]
    }
    
    # Initialize annotator
    annotator = ClinVarAnnotator()
    
    # Annotate the variant
    result = annotator.annotate_variant(test_variant)
    
    # Print results
    print(f"\n{'='*80}")
    print("CLINVAR ANNOTATION RESULTS")
    print(f"{'='*80}")
    print(f"Variant ID: {result.variant_id}")
    print(f"Gene: {result.gene}")
    print(f"Protein Change: {result.protein_change}")
    print(f"Clinical Significance: {result.clinical_significance}")
    print(f"Review Status: {result.review_status}")
    print(f"Condition: {result.condition}")
    print(f"Variation ID: {result.variation_id}")
    print(f"Last Evaluated: {result.last_evaluated}")
    print(f"Submitter: {result.submitter}")
    print(f"Star Rating: {result.star_rating}")
    print(f"Accession: {result.accession}")
    print(f"Source: {result.source}")
    
    if result.warnings:
        print(f"\nWarnings:")
        for warning in result.warnings:
            print(f"  - {warning}")
    
    print(f"\n{'='*80}")


if __name__ == "__main__":
    main() 