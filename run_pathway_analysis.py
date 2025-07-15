#!/usr/bin/env python3
"""
Command-line script to run pathway impact analysis
"""

import sys
import argparse
from pathlib import Path

# Add modules to path
sys.path.append(str(Path(__file__).parent / "modules"))

from pathway_impact_analyzer import PathwayImpactAnalyzer


def main():
    """Main function with command line arguments"""
    parser = argparse.ArgumentParser(
        description="Run pathway impact analysis using GSEApy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pathway_analysis.py
  python run_pathway_analysis.py -i variant_impact_results_large.json
  python run_pathway_analysis.py -i my_variants.json -o my_results.json
        """
    )
    
    parser.add_argument(
        "-i", "--input", 
        default="variant_impact_results.json",
        help="Input JSON file with variant data (default: variant_impact_results.json)"
    )
    
    parser.add_argument(
        "-o", "--output", 
        default="pathway_impact_results.json",
        help="Output JSON file for results (default: pathway_impact_results.json)"
    )
    
    parser.add_argument(
        "-m", "--method",
        default="pathogenicity_mean",
        choices=["pathogenicity_mean", "pathogenicity_max", "pathogenicity_sum", "uniform", "variant_count", "clinical_relevance"],
        help="Gene ranking method (default: pathogenicity_mean)"
    )
    
    parser.add_argument(
        "--gene-sets",
        nargs="+",
        default=["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2021"],
        help="Gene set databases to use (default: KEGG_2021_Human Reactome_2022 GO_Biological_Process_2021)"
    )
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.input).exists():
        print(f"❌ Error: Input file '{args.input}' not found")
        return 1
    
    print(f"🔬 Running Pathway Impact Analysis")
    print(f"   Input: {args.input}")
    print(f"   Output: {args.output}")
    print(f"   Method: {args.method}")
    print(f"   Gene sets: {', '.join(args.gene_sets)}")
    print("=" * 60)
    
    try:
        # Initialize analyzer
        analyzer = PathwayImpactAnalyzer(input_file=args.input)
        
        # Load data and extract genes
        if not analyzer.load_input_data():
            print("❌ Failed to load input data")
            return 1
        
        genes = analyzer.extract_gene_list()
        if not genes:
            print("❌ No genes found in input data")
            return 1
        
        print(f"📊 Found {len(genes)} genes: {', '.join(genes[:10])}{'...' if len(genes) > 10 else ''}")
        
        # Create gene ranking
        gene_scores = analyzer.create_gene_ranking(method=args.method)
        
        # Run GSEA analysis
        enrichment_results = analyzer.run_gsea_analysis(gene_scores, gene_sets=args.gene_sets)
        
        # Analyze target pathways
        pathway_scores = analyzer.analyze_target_pathways()
        
        # Generate report
        report = analyzer.generate_report(output_file=args.output)
        
        # Print results
        print("\n" + "="*60)
        print("🎯 TARGET PATHWAY VALIDATION")
        print("="*60)
        
        validation_results = analyzer.validate_target_scores()
        all_passed = True
        
        for pathway, validation in validation_results.items():
            status = validation['status']
            score = validation['actual_score']
            expected = validation['expected_min']
            passed = validation['passed']
            
            if not passed:
                all_passed = False
            
            icon = "✅" if passed else "❌"
            print(f"{icon} {pathway:15} | {score:8.3f} | {'≥' if passed else '<':2} {expected:6.2f} | {status}")
        
        print("\n" + "="*60)
        if all_passed:
            print("🎉 ALL TARGET PATHWAYS PASSED VALIDATION!")
        else:
            print("⚠️  SOME TARGET PATHWAYS FAILED VALIDATION")
        print("="*60)
        
        print(f"\n📄 Results saved to: {args.output}")
        print("✅ Analysis completed successfully!")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main()) 