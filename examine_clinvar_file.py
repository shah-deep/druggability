#!/usr/bin/env python3
"""
Examine ClinVar variant summary file structure
"""

import gzip
import sys

def examine_clinvar_file():
    """Examine the structure of the ClinVar variant summary file"""
    
    file_path = "cache/clinvar_variant_summary.txt.gz"
    
    try:
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            # Read header
            header = f.readline().strip()
            print("Header:")
            print(header)
            print(f"\nHeader fields: {header.split('\t')}")
            print(f"Number of fields: {len(header.split('\t'))}")
            
            # Read first few data lines
            print("\nFirst 3 data lines:")
            for i in range(3):
                line = f.readline().strip()
                if line:
                    fields = line.split('\t')
                    print(f"Line {i+1}: {len(fields)} fields")
                    print(f"  Gene: {fields[0] if len(fields) > 0 else 'N/A'}")
                    print(f"  Name: {fields[1] if len(fields) > 1 else 'N/A'}")
                    print(f"  ClinicalSignificance: {fields[2] if len(fields) > 2 else 'N/A'}")
                    print(f"  VariationID: {fields[3] if len(fields) > 3 else 'N/A'}")
                    print()
            
            # Count total lines
            print("Counting total lines...")
            f.seek(0)  # Reset to beginning
            next(f)  # Skip header
            line_count = sum(1 for line in f)
            print(f"Total data lines: {line_count:,}")
            
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Look for specific variants
    print("\nSearching for specific variants...")
    search_variants = [
        ("TP53", "p.Arg175His"),
        ("TP53", "p.Arg273His"),
        ("CFTR", "p.Phe508del"),
        ("BRCA1", "p.Cys61Gly")
    ]
    
    try:
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            header = f.readline().strip()
            header_fields = header.split('\t')
            
            # Find column indices
            gene_idx = header_fields.index('GeneSymbol') if 'GeneSymbol' in header_fields else 0
            name_idx = header_fields.index('Name') if 'Name' in header_fields else 1
            sig_idx = header_fields.index('ClinicalSignificance') if 'ClinicalSignificance' in header_fields else 2
            var_id_idx = header_fields.index('VariationID') if 'VariationID' in header_fields else 3
            
            print(f"Column indices - Gene: {gene_idx}, Name: {name_idx}, Significance: {sig_idx}, VariationID: {var_id_idx}")
            
            for gene, protein_change in search_variants:
                print(f"\nSearching for {gene} {protein_change}...")
                f.seek(0)
                next(f)  # Skip header
                
                found = False
                for line_num, line in enumerate(f, 1):
                    fields = line.strip().split('\t')
                    if len(fields) > max(gene_idx, name_idx, sig_idx, var_id_idx):
                        file_gene = fields[gene_idx]
                        file_name = fields[name_idx]
                        
                        if file_gene == gene and protein_change in file_name:
                            significance = fields[sig_idx]
                            variation_id = fields[var_id_idx]
                            print(f"  Found: {gene} {protein_change}")
                            print(f"    Clinical Significance: {significance}")
                            print(f"    Variation ID: {variation_id}")
                            print(f"    Full Name: {file_name}")
                            found = True
                            break
                
                if not found:
                    print(f"  Not found: {gene} {protein_change}")
                    
    except Exception as e:
        print(f"Error searching file: {e}")

if __name__ == "__main__":
    examine_clinvar_file() 