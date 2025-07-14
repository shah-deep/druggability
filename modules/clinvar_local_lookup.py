#!/usr/bin/env python3
"""
ClinVar Local Lookup

Indexes and looks up clinical significance from the ClinVar variant_summary.txt file.
"""
import gzip
import os
from typing import Dict, Optional, List

class ClinVarLocalLookup:
    def __init__(self, summary_path: str = 'cache/clinvar_variant_summary.txt.gz'):
        self.summary_path = summary_path
        self.index = {}
        if os.path.exists(self.summary_path):
            self._build_index()

    def _build_index(self):
        with gzip.open(self.summary_path, 'rt', encoding='utf-8') as f:
            header = next(f).strip().split('\t')
            col_idx = {col: i for i, col in enumerate(header)}
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < len(header):
                    continue
                gene = fields[col_idx['GeneSymbol']]
                protein = fields[col_idx['Name']]
                significance = fields[col_idx['ClinicalSignificance']]
                variation_id = fields[col_idx['VariationID']]
                # Index by gene + protein_change (HGVS or protein string)
                key = f"{gene}|{protein}"
                self.index[key] = {
                    'significance': significance,
                    'variation_id': variation_id,
                    'raw': fields
                }

    def lookup(self, gene: str, protein_change: str) -> Optional[Dict]:
        # Try exact match
        key = f"{gene}|{protein_change}"
        if key in self.index:
            return self.index[key]
        # Try partial match (protein_change substring)
        for k, v in self.index.items():
            if k.startswith(f"{gene}|") and protein_change in k:
                return v
        return None

    def lookup_by_variation_id(self, variation_id: str) -> Optional[Dict]:
        for v in self.index.values():
            if v['variation_id'] == variation_id:
                return v
        return None 