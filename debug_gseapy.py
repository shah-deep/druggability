#!/usr/bin/env python3
"""
Debug script to test GSEApy functionality
"""

import sys
import os
import json
import logging
from pathlib import Path

# Add modules to path
sys.path.append(str(Path(__file__).parent / "modules"))

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def test_gseapy_installation():
    """Test if GSEApy is working"""
    try:
        import gseapy as gp
        from gseapy import GSEA, enrichr
        print("✅ GSEApy imported successfully")
        
        # Test available gene sets
        print("\nTesting available gene sets...")
        try:
            # Try to get available gene sets
            gs = gp.get_library_name()
            print(f"Available gene set libraries: {gs[:10]}...")  # Show first 10
            return True
        except Exception as e:
            print(f"❌ Error getting gene sets: {e}")
            return False
            
    except ImportError as e:
        print(f"❌ GSEApy import failed: {e}")
        return False


def test_enrichr_with_simple_data():
    """Test enrichr with simple gene list"""
    try:
        from gseapy import enrichr
        
        # Simple test gene list
        test_genes = ['TP53', 'BRCA1', 'BRCA2', 'ATM', 'CHEK2']
        
        print(f"\nTesting enrichr with genes: {test_genes}")
        
        # Try with KEGG
        enr = enrichr(
            gene_list=test_genes,
            gene_sets=['KEGG_2021_Human'],
            organism='Human',
            outdir=None,
            no_plot=True
        )
        
        if hasattr(enr, 'results') and enr.results is not None:
            print(f"✅ Enrichr worked! Found {len(enr.results)} results")
            print("Sample results:")
            print(enr.results.head(3))
            return True
        else:
            print("❌ Enrichr returned no results")
            return False
            
    except Exception as e:
        print(f"❌ Enrichr test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_with_actual_data():
    """Test with actual data from variant file"""
    try:
        # Load variant data
        with open("variant_impact_results_large.json", 'r') as f:
            data = json.load(f)
        
        # Extract genes
        genes = set()
        for variant in data.get('missense_variants', []):
            if 'gene' in variant and variant['gene']:
                genes.add(variant['gene'])
        
        gene_list = list(genes)
        print(f"\nTesting with actual genes: {gene_list[:10]}... (total: {len(gene_list)})")
        
        from gseapy import enrichr
        
        # Try with KEGG
        enr = enrichr(
            gene_list=gene_list,
            gene_sets=['KEGG_2021_Human'],
            organism='Human',
            outdir=None,
            no_plot=True
        )
        
        if hasattr(enr, 'results') and enr.results is not None and not enr.results.empty:
            print(f"✅ Enrichr worked with actual data! Found {len(enr.results)} results")
            print("Sample results:")
            print(enr.results.head(5))
            
            # Check for target pathways
            target_pathways = ['p53', 'DNA repair', 'apoptosis']
            for pathway in target_pathways:
                matches = enr.results[enr.results['Term'].str.contains(pathway, case=False, na=False)]
                if not matches.empty:
                    print(f"\nFound pathways containing '{pathway}':")
                    print(matches[['Term', 'P-value', 'Adjusted P-value']].head(3))
                else:
                    print(f"\nNo pathways found containing '{pathway}'")
            
            return True
        else:
            print("❌ Enrichr returned no results with actual data")
            return False
            
    except Exception as e:
        print(f"❌ Test with actual data failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main debug function"""
    print("🔍 Debugging GSEApy Analysis")
    print("=" * 50)
    
    # Test 1: Installation
    print("\n1. Testing GSEApy installation...")
    if not test_gseapy_installation():
        return
    
    # Test 2: Simple enrichr test
    print("\n2. Testing enrichr with simple data...")
    if not test_enrichr_with_simple_data():
        return
    
    # Test 3: Test with actual data
    print("\n3. Testing with actual variant data...")
    test_with_actual_data()


if __name__ == "__main__":
    main() 