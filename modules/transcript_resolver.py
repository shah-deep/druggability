#!/usr/bin/env python3
"""
Enhanced Transcript ID Resolver with VEP API Integration

This module provides robust functionality to resolve exact transcript IDs for specific variants.
Uses Ensembl VEP API with multiple fallback strategies for comprehensive transcript resolution.
"""

import requests
import json
import logging
import time
from typing import Dict, List, Optional, Tuple
from pathlib import Path
from urllib.parse import quote

# Configure logging
from .logging_config import setup_logging, get_logger
setup_logging()
logger = get_logger(__name__)

class TranscriptResolver:
    """Resolve exact transcript IDs for specific variants using Ensembl VEP API"""
    
    def __init__(self, cache_dir: str = "cache"):
        self.vep_api_base = "https://rest.ensembl.org/vep"
        self.ensembl_api_base = "https://rest.ensembl.org"
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.cache_file = self.cache_dir / "transcript_mappings.json"
        self.cache = {} #self._load_cache()
        
        # API rate limiting
        self.request_delay = 0.05  # 100ms between requests
        self.max_retries = 1
        self.timeout = 30
        
        # Headers for API requests
        self.headers = {
            "Content-Type": "application/json",
            "User-Agent": "TranscriptResolver/1.0"
        }
    
    def _load_cache(self) -> Dict:
        """Load transcript mapping cache"""
        try:
            if self.cache_file.exists():
                with open(self.cache_file, 'r') as f:
                    return json.load(f)
        except Exception as e:
            logger.warning(f"Could not load transcript cache: {e}")
        return {}
    
    def _save_cache(self):
        """Save transcript mapping cache"""
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save transcript cache: {e}")
    
    def _make_api_request(self, url: str, method: str = "GET", data: Optional[Dict] = None) -> Optional[Dict]:
        """Make API request with rate limiting and retry logic"""
        for attempt in range(self.max_retries):
            try:
                time.sleep(self.request_delay)  # Rate limiting
                
                if method == "GET":
                    response = requests.get(url, headers=self.headers, timeout=self.timeout)
                elif method == "POST":
                    response = requests.post(url, headers=self.headers, json=data, timeout=self.timeout)
                else:
                    logger.error(f"Unsupported HTTP method: {method}")
                    return None
                
                if response.status_code == 200:
                    return response.json()
                elif response.status_code == 429:  # Rate limited
                    wait_time = (attempt + 1) * 2
                    logger.warning(f"Rate limited, waiting {wait_time}s before retry")
                    time.sleep(wait_time)
                    continue
                else:
                    logger.warning(f"API request failed: {response.status_code} - {response.text}")
                    return None
                    
            except requests.exceptions.RequestException as e:
                logger.warning(f"Request failed (attempt {attempt + 1}): {e}")
                if attempt == self.max_retries - 1:
                    return None
                time.sleep(1)  # Wait before retry
        
        return None

    def _prepare_variant_for_vep(self, variant: Dict) -> Optional[Dict]:
        """Fill missing chromosome, normalize alleles for VEP, return new dict or None if not possible."""
        v = variant.copy()
        
        # 1. Chromosome resolution
        chrom = v.get('chromosome')
        if not chrom:
            gene = v.get('gene')
            if gene:
                chrom = self._fetch_chromosome_from_gene(gene)
                if chrom:
                    v['chromosome'] = chrom
                    logger.info(f"Resolved chromosome for {gene}: {chrom}")
                else:
                    logger.warning(f"Could not determine chromosome for gene {gene}")
                    return None
            else:
                logger.warning(f"Variant missing chromosome and gene: {variant}")
                return None
        
        # 2. Normalize alleles for VEP
        ref = v.get('reference', '')
        alt = v.get('alternate', '')
        
        # Handle deletions (empty alternate)
        if ref and not alt:
            v['alternate'] = '-'
            logger.info(f"Normalized deletion: {ref} -> -")
        
        # Handle insertions (empty reference)
        if alt and not ref:
            v['reference'] = '-'
            logger.info(f"Normalized insertion: - -> {alt}")
        
        # 3. Validate position
        if not v.get('position', 0):
            logger.warning(f"Missing position for variant: {variant}")
            return None
        
        return v

    def _fetch_chromosome_from_gene(self, gene: str) -> Optional[str]:
        """Fetch chromosome for a gene using Ensembl API."""
        try:
            url = f"{self.ensembl_api_base}/lookup/symbol/homo_sapiens/{quote(gene)}"
            data = self._make_api_request(url)
            
            if data and 'seq_region_name' in data:
                return str(data['seq_region_name'])
            
        except Exception as e:
            logger.warning(f"Could not fetch chromosome for {gene}: {e}")
        return None

    def resolve_transcript_id(self, gene: str, protein_change: str, 
                            chromosome: str = "", position: int = 0,
                            reference: str = "", alternate: str = "") -> Optional[str]:
        """Resolve exact transcript ID for a specific variant"""
        
        # Prepare variant dict
        variant = {
            'gene': gene,
            'protein_change': protein_change,
            'chromosome': chromosome,
            'position': position,
            'reference': reference,
            'alternate': alternate
        }
        
        v = self._prepare_variant_for_vep(variant)
        if not v:
            logger.warning(f"Could not prepare variant for VEP: {variant}")
            return None
        
        # Create cache key
        cache_key = f"{v['gene']}_{v['protein_change']}_{v['chromosome']}_{v['position']}_{v['reference']}_{v['alternate']}"
        
        # Check cache first
        if cache_key in self.cache:
            logger.info(f"Cache hit for {gene} {protein_change}")
            return self.cache[cache_key]
        
        # Try multiple VEP strategies
        transcript_id = self._resolve_with_vep_strategies(v)
        
        if transcript_id:
            self.cache[cache_key] = transcript_id
            # self._save_cache()
            logger.info(f"Resolved transcript ID for {gene} {protein_change}: {transcript_id}")
        else:
            logger.warning(f"Could not resolve transcript ID for {gene} {protein_change}")
        
        return transcript_id

    def get_gencode_id(self, gene: str) -> Optional[str]:
        """Get the gencode_id for a gene by constructing id.version from transcript data"""
        try:
            url = f"{self.ensembl_api_base}/lookup/symbol/homo_sapiens/{quote(gene)}?expand=1"
            data = self._make_api_request(url)
            
            if not data or 'id' not in data:
                return None
            
            id = data['id']
            version = data['version']
            gencode_id = f"{id}.{version}"
            logger.info(f"Found gencode_id for {gene}: {gencode_id}")
            return gencode_id
            
        except Exception as e:
            logger.warning(f"Error getting gencode_id for {gene}: {e}")
            return None
    
    def _resolve_with_vep_strategies(self, variant: Dict) -> Optional[str]:
        """Try multiple VEP strategies to resolve transcript ID"""
        
        gene = variant['gene']
        
        # Strategy 1: VEP with genomic coordinates (most reliable)
        # transcript_id = self._query_vep_for_variant(variant)
        # if transcript_id:
        #     return transcript_id
        
        # Strategy 2: Gene lookup for canonical transcript (fallback)
        logger.info(f"Falling back to canonical transcript for {gene}")
        transcript_id = self._get_canonical_transcript_from_gene(gene)
        if transcript_id:
            return transcript_id

        return None

    def _query_vep_for_variant(self, variant: Dict) -> Optional[str]:
        """Query VEP with a variant dictionary using gene symbol and protein change"""
        gene = variant.get('gene')
        protein_change = variant.get('protein_change')
        chromosome = variant.get('chromosome')
        position = variant.get('position')
        reference = variant.get('reference')
        alternate = variant.get('alternate')

        if not gene or not protein_change:
            logger.warning(f"Missing gene or protein_change for variant: {variant}")
            return None

        # Strategy 1: Use gene symbol lookup first (most reliable for known genes)
        transcript_id = self._query_vep_by_gene_symbol(gene, protein_change)
        if transcript_id:
            return transcript_id
        
        # Strategy 2: Fall back to genomic coordinates if available
        if all([chromosome, position, reference, alternate]):
            # Ensure proper types
            chrom_str = str(chromosome) if chromosome else ""
            pos_int = int(position) if position else 0
            ref_str = str(reference) if reference else ""
            alt_str = str(alternate) if alternate else ""
            
            if chrom_str and pos_int and ref_str and alt_str:
                transcript_id = self._query_vep_by_genomic_coordinates(chrom_str, pos_int, ref_str, alt_str, gene, protein_change)
                if transcript_id:
                    return transcript_id
        
        return None

    def _query_vep_by_gene_symbol(self, gene: str, protein_change: str) -> Optional[str]:
        """Query VEP using gene symbol and protein change"""
        try:
            logger.info(f"Querying VEP with gene symbol: {gene} and protein change: {protein_change}")
            
            # Use Ensembl gene lookup to get gene information
            url = f"{self.ensembl_api_base}/lookup/symbol/homo_sapiens/{quote(gene)}?expand=1"
            gene_data = self._make_api_request(url)
            
            if not gene_data or 'Transcript' not in gene_data:
                logger.warning(f"No gene data found for {gene}")
                return None
            
            # Look for transcript with matching protein change
            for transcript in gene_data['Transcript']:
                transcript_id = transcript.get('id')
                if transcript_id:
                    # Check if this transcript has the protein change
                    if self._check_transcript_protein_change(transcript_id, protein_change):
                        logger.info(f"Found transcript {transcript_id} with protein change {protein_change}")
                        return transcript_id
            
            # If no exact match, return canonical transcript
            for transcript in gene_data['Transcript']:
                if transcript.get('is_canonical', 0) == 1:
                    transcript_id = transcript.get('id')
                    logger.info(f"Using canonical transcript for {gene}: {transcript_id}")
                    return transcript_id
            
            # Fallback to first transcript
            if gene_data['Transcript']:
                transcript_id = gene_data['Transcript'][0].get('id')
                logger.info(f"Using first transcript for {gene}: {transcript_id}")
                return transcript_id
            
        except Exception as e:
            logger.warning(f"Error querying VEP by gene symbol: {e}")
        
        return None

    def _check_transcript_protein_change(self, transcript_id: str, protein_change: str) -> bool:
        """Check if a transcript has the specified protein change using VEP"""
        try:
            # Query VEP for this specific transcript
            url = f"{self.vep_api_base}/homo_sapiens/id/{transcript_id}"
            data = self._make_api_request(url)
            
            if data and isinstance(data, list):
                for entry in data:
                    if 'transcript_consequences' in entry:
                        for transcript in entry['transcript_consequences']:
                            if transcript.get('transcript_id') == transcript_id:
                                hgvsp = transcript.get('hgvsp', '')
                                if hgvsp and protein_change in hgvsp:
                                    return True
                                # Also check amino_acids field
                                amino_acids = transcript.get('amino_acids', '')
                                if amino_acids and protein_change in amino_acids:
                                    return True
            
        except Exception as e:
            logger.warning(f"Error checking protein change for transcript {transcript_id}: {e}")
        
        return False

    def _query_vep_by_genomic_coordinates(self, chromosome: str, position: int, 
                                         reference: str, alternate: str, 
                                         gene: str, protein_change: str) -> Optional[str]:
        """Query VEP API using genomic coordinates as fallback"""
        try:
            variant_str = f"{chromosome} {position} . {reference} {alternate} . . ."
            logger.info(f"Querying VEP with genomic coordinates: {variant_str}")
            
            url = f"{self.vep_api_base}/homo_sapiens/region"
            payload = {"variants": [variant_str]}
            data = self._make_api_request(url, method="POST", data=payload)
            
            if data and isinstance(data, list) and data:
                return self._extract_transcript_from_vep(data, gene, protein_change)
        
        except Exception as e:
            logger.warning(f"Error querying VEP with genomic coordinates: {e}")
        
        return None

    def _construct_hgvs_notation(self, chromosome: str, position: int, reference: str, alternate: str) -> Optional[str]:
        """Construct HGVS notation from genomic coordinates"""
        try:
            # Basic HGVS format: NC_0000XX.X:g.POSITIONREF>ALT
            # For now, return None as this requires more complex logic
            # This could be enhanced with proper chromosome mapping
            return None
        except Exception as e:
            logger.warning(f"Could not construct HGVS notation: {e}")
            return None
    
    def _get_canonical_transcript_from_gene(self, gene: str) -> Optional[str]:
        """Get canonical transcript ID using Ensembl gene lookup"""
        try:
            url = f"{self.ensembl_api_base}/lookup/symbol/homo_sapiens/{quote(gene)}?expand=1"
            data = self._make_api_request(url)
            
            if not data:
                return None
            
            # First, try to get the canonical transcript from the gene-level canonical_transcript field
            canonical_transcript = data.get('canonical_transcript')
            if canonical_transcript:
                logger.info(f"Found canonical transcript for {gene}: {canonical_transcript}")
                return canonical_transcript
            
            # Fallback: Look through transcripts for canonical
            if 'Transcript' in data:
                for transcript in data['Transcript']:
                    if transcript.get('is_canonical', 0) == 1:
                        transcript_id = transcript.get('id')
                        version = transcript.get('version')
                        # Return stable ID with decimal format if version is available
                        if version is not None:
                            stable_id_with_version = f"{transcript_id}.{version}"
                            logger.info(f"Found canonical transcript for {gene}: {stable_id_with_version}")
                            return stable_id_with_version
                        else:
                            logger.info(f"Found canonical transcript for {gene}: {transcript_id}")
                            return transcript_id
                
                # If no canonical, return first transcript
                if data['Transcript']:
                    transcript_id = data['Transcript'][0].get('id')
                    version = data['Transcript'][0].get('version')
                    # Return stable ID with decimal format if version is available
                    if version is not None:
                        stable_id_with_version = f"{transcript_id}.{version}"
                        logger.info(f"Found first transcript for {gene}: {stable_id_with_version}")
                        return stable_id_with_version
                    else:
                        logger.info(f"Found first transcript for {gene}: {transcript_id}")
                        return transcript_id
            
            return None
            
        except Exception as e:
            logger.warning(f"Error in gene lookup: {e}")
            return None
    
    def _extract_transcript_from_vep(self, vep_data: List[Dict], gene: str, protein_change: Optional[str]) -> Optional[str]:
        """Extract exact transcript ID from VEP response"""
        
        all_consequences = []
        for variant_data in vep_data:
            if 'transcript_consequences' in variant_data:
                all_consequences.extend(variant_data['transcript_consequences'])

        if not all_consequences:
            logger.warning(f"No transcript_consequences found for {gene} {protein_change}")
            return None

        # Sort transcripts by preference
        sorted_transcripts = self._sort_transcripts(all_consequences)

        # Find the best matching transcript
        for transcript in sorted_transcripts:
            if self._is_correct_transcript(transcript, gene, protein_change):
                transcript_id = transcript.get('transcript_id')
                logger.info(f"Found matching transcript for {gene} {protein_change}: {transcript_id}")
                return transcript_id
        
        logger.warning(f"No suitable transcript found for {gene} {protein_change} after filtering.")
        # Fallback to the top-sorted transcript if no perfect match is found
        if sorted_transcripts:
            top_transcript_id = sorted_transcripts[0].get('transcript_id')
            logger.info(f"Falling back to top-sorted transcript: {top_transcript_id}")
            return top_transcript_id

        return None

    def _sort_transcripts(self, transcripts: List[Dict]) -> List[Dict]:
        """Sort transcripts by preference: MANE, Canonical, APPRIS, TSL."""
        def get_score(t):
            score = 0
            if t.get('mane_select'): score += 10
            if t.get('canonical'): score += 5
            if t.get('appris'):
                if 'principal' in t['appris']: score += 2
            # Lower TSL is better (Transcript Support Level)
            tsl = t.get('tsl')
            if tsl is not None:
                score += (5 - tsl) # tsl 1 is best (score 4), tsl 5 is worst (score 0)
            return score

        return sorted(transcripts, key=get_score, reverse=True)

    def _is_correct_transcript(self, transcript: Dict, gene: str, protein_change: Optional[str]) -> bool:
        """Check if transcript is correct for the specific variant."""
        
        # 1. Check gene name
        transcript_gene = transcript.get('gene_symbol', '')
        if transcript_gene.upper() != gene.upper():
            return False
        
        # 2. Check protein change if available (relaxed check)
        if protein_change:
            hgvsp = transcript.get('hgvsp', '')
            # e.g., protein_change "p.Arg175His" should be in hgvsp "ENST...:p.Arg175His"
            if hgvsp and protein_change in hgvsp:
                return True
        
        # 3. If no matching protein change, check consequence terms as a fallback
        consequence_terms = transcript.get('consequence_terms', [])
        valid_consequences = [
            'missense_variant', 
            'inframe_deletion', 
            'frameshift_variant', 
            'inframe_insertion'
        ]
        if any(term in valid_consequences for term in consequence_terms):
            return True

        return False

    def _is_correct_protein_consequence(self, protein: Dict, gene: str, protein_change: str) -> bool:
        """Check if protein consequence is correct for the specific variant"""
        
        # Check gene name
        protein_gene = protein.get('gene_symbol', '')
        if protein_gene.upper() != gene.upper():
            return False
        
        # Check protein change if available
        if 'hgvsp' in protein and protein['hgvsp'].endswith(protein_change):
            return True

        if 'protein_change' in protein and protein['protein_change'] == protein_change:
            return True
        
        return False
    
    def resolve_variants(self, variants: List[Dict]) -> List[Dict]:
        """Resolve transcript IDs for a list of variants and include VEP transcript IDs"""
        resolved_variants = []
        
        for variant in variants:
            v = self._prepare_variant_for_vep(variant)
            if not v:
                resolved_variants.append(variant)
                continue
            
            # Extract variant information
            gene = v.get('gene', '')
            protein_change = v.get('protein_change', '')
            chromosome = v.get('chromosome', '')
            position = v.get('position', 0)
            reference = v.get('reference', '')
            alternate = v.get('alternate', '')
            
            # Resolve transcript ID
            transcript_id = None
            if gene and protein_change and not variant.get('transcript_id'):
                transcript_id = self.resolve_transcript_id(
                    gene, protein_change, chromosome, position, reference, alternate
                )
                if transcript_id:
                    v['transcript_id'] = transcript_id
            
            # Always fetch VEP transcript IDs for this variant (using POST endpoint)
            vep_transcript_ids = self._get_all_vep_transcripts(v)
            v['vep_transcript_ids'] = vep_transcript_ids
            
            resolved_variants.append(v)
        
        return resolved_variants

    def _get_all_vep_transcripts(self, variant: Dict) -> List[str]:
        """Fetch all transcript IDs from VEP for a given variant using gene symbol"""
        try:
            gene = variant.get('gene')
            protein_change = variant.get('protein_change')
            
            if not gene:
                return []

            # Use gene symbol lookup to get all transcripts
            url = f"{self.ensembl_api_base}/lookup/symbol/homo_sapiens/{quote(gene)}?expand=1"
            gene_data = self._make_api_request(url)
            
            if not gene_data or 'Transcript' not in gene_data:
                return []

            transcript_ids = []
            for transcript in gene_data['Transcript']:
                transcript_id = transcript.get('id')
                if transcript_id:
                    transcript_ids.append(transcript_id)
            
            logger.info(f"Found {len(transcript_ids)} transcripts for gene {gene}")
            return transcript_ids
            
        except Exception as e:
            logger.warning(f"Error fetching VEP transcript IDs for gene {gene}: {e}")
            return []

    def get_cache_stats(self) -> Dict:
        """Get cache statistics"""
        return {
            'cache_size': len(self.cache),
            'cache_file': str(self.cache_file),
            'cache_dir': str(self.cache_dir)
        }

    def clear_cache(self):
        """Clear the transcript cache"""
        self.cache = {}
        if self.cache_file.exists():
            self.cache_file.unlink()
        logger.info("Transcript cache cleared") 