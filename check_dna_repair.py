#!/usr/bin/env python3
"""
Check for DNA repair pathways in enrichment results
"""

import json

def check_dna_repair():
    with open('pathway_impact_results.json', 'r') as f:
        data = json.load(f)
    
    results = data['enrichment_results']
    dna_repair_paths = []
    
    for gene_set, gs_data in results.items():
        if 'enrichment_results' in gs_data:
            for pathway in gs_data['enrichment_results']:
                term = pathway.get('Term', '').lower()
                if 'dna' in term or 'repair' in term:
                    dna_repair_paths.append({
                        'gene_set': gene_set,
                        'term': pathway.get('Term', ''),
                        'es': pathway.get('ES', 0),
                        'nes': pathway.get('NES', 0),
                        'fdr': pathway.get('FDR q-val', 1.0)
                    })
    
    print(f"Found {len(dna_repair_paths)} DNA repair related pathways:")
    for i, pathway in enumerate(dna_repair_paths):
        print(f"  {i+1}. {pathway['term']}")
        print(f"     Gene set: {pathway['gene_set']}")
        print(f"     ES: {pathway['es']:.3f}, NES: {pathway['nes']:.3f}, FDR: {pathway['fdr']:.3e}")
        print()
    
    if not dna_repair_paths:
        print("No DNA repair pathways found!")
        print("\nChecking for similar terms...")
        
        # Look for similar terms
        similar_terms = []
        for gene_set, gs_data in results.items():
            if 'enrichment_results' in gs_data:
                for pathway in gs_data['enrichment_results']:
                    term = pathway.get('Term', '').lower()
                    if any(word in term for word in ['homologous', 'recombination', 'mismatch', 'nucleotide', 'base excision']):
                        similar_terms.append({
                            'gene_set': gene_set,
                            'term': pathway.get('Term', ''),
                            'es': pathway.get('ES', 0),
                            'nes': pathway.get('NES', 0),
                            'fdr': pathway.get('FDR q-val', 1.0)
                        })
        
        print(f"Found {len(similar_terms)} related pathways:")
        for i, pathway in enumerate(similar_terms[:5]):
            print(f"  {i+1}. {pathway['term']}")
            print(f"     Gene set: {pathway['gene_set']}")
            print(f"     ES: {pathway['es']:.3f}, NES: {pathway['nes']:.3f}, FDR: {pathway['fdr']:.3e}")
            print()

if __name__ == "__main__":
    check_dna_repair() 