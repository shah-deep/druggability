# CADD Pathogenicity Scoring Guide

## Overview

This guide explains how to use the **CADD (Combined Annotation Dependent Depletion)** pathogenicity scorer, which is the **gold standard** method for variant pathogenicity assessment in clinical genomics.

## What is CADD?

CADD is a widely accepted and standardized method for predicting the pathogenicity of genetic variants. It integrates multiple annotations including:

- **Conservation scores** across species
- **Regulatory annotations** (promoters, enhancers, etc.)
- **Protein function predictions**
- **Disease associations** from databases
- **Population frequencies** from large cohorts

## Why CADD is the Gold Standard

1. **Comprehensive Integration**: Combines multiple data sources into a single score
2. **Clinically Validated**: Widely used in clinical interpretation pipelines
3. **Standardized Scoring**: CADD Phred scores provide interpretable pathogenicity levels
4. **Regular Updates**: Continuously updated with new annotations and data

## CADD Score Interpretation

### CADD Phred Score Ranges

| Score Range | Pathogenicity Level | Clinical Interpretation |
|-------------|-------------------|------------------------|
| ≥ 30 | HIGH | Likely deleterious |
| 20-30 | MEDIUM | Possibly deleterious |
| 10-20 | LOW | Possibly tolerated |
| < 10 | BENIGN | Likely tolerated |

### CADD Raw Score
- Range: 0 to 100
- Higher values indicate more deleterious variants
- Used internally by CADD algorithm

## Installation and Setup

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Create Required Directories

```bash
mkdir -p cache/cadd
mkdir -p data/cadd
```

### 3. Optional: Download CADD Data

For offline analysis, you can download pre-computed CADD scores:

```bash
# Download CADD scores for hg38 (recommended)
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz -O data/cadd/cadd_scores.tsv.gz

# Or download for hg19
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz -O data/cadd/cadd_scores.tsv.gz
```

## Usage Examples

### Basic CADD Scoring

```python
from modules.cadd_pathogenicity_scorer import PathogenicityScorer

# Initialize scorer
scorer = PathogenicityScorer()

# Define variants
variants = [
    {
        'id': 'rs121913343',
        'chromosome': '7',
        'position': 140453136,
        'reference_allele': 'A',
        'alternate_allele': 'T',
        'gene': 'BRAF',
        'protein_change': 'p.Val600Glu'
    }
]

# Score variants
results = scorer.score_variants(variants)

# Access results
for variant_id, score in results.items():
    print(f"Variant: {variant_id}")
    print(f"CADD Phred Score: {score.cadd_phred_score:.2f}")
    print(f"Pathogenicity Level: {score.pathogenicity_level}")
    print(f"Is Pathogenic: {score.is_pathogenic}")
```

### Comprehensive Analysis

```python
from modules.integrated_pathogenicity_analyzer import IntegratedPathogenicityAnalyzer

# Initialize integrated analyzer
analyzer = IntegratedPathogenicityAnalyzer()

# Perform comprehensive analysis
results = analyzer.analyze_variants_comprehensive(variants)

# Export results
analyzer.export_comprehensive_results(results, 'pathogenicity_analysis.json')
```

### Command Line Usage

```bash
# Basic CADD scoring
python example_cadd_usage.py

# Comprehensive analysis
python modules/integrated_pathogenicity_analyzer.py
```

## Input Format

### Required Fields

```python
variant = {
    'id': 'rs121913343',                    # Variant identifier
    'chromosome': '7',                      # Chromosome number
    'position': 140453136,                  # Genomic position
    'reference_allele': 'A',                # Reference allele
    'alternate_allele': 'T',                # Alternate allele
    'gene': 'BRAF',                         # Gene name (optional)
    'protein_change': 'p.Val600Glu'         # Protein change (optional)
}
```

### Supported Chromosome Formats

- `'7'` or `'chr7'` → `'7'`
- `'X'` or `'chrX'` → `'X'`
- `'Y'` or `'chrY'` → `'Y'`

## Output Format

### CADD Score Object

```python
CADDScore(
    variant_id='rs121913343',
    cadd_raw_score=4.123,
    cadd_phred_score=25.67,
    chromosome='7',
    position=140453136,
    reference_allele='A',
    alternate_allele='T',
    gene='BRAF',
    transcript='ENST00000288602',
    protein_change='p.Val600Glu',
    consequence='missense_variant',
    pathogenicity_level='MEDIUM',
    is_pathogenic=True,
    prediction_source='CADD API'
)
```

### JSON Export Format

```json
{
  "metadata": {
    "timestamp": "2025-01-27T10:30:00",
    "total_variants": 1,
    "scoring_method": "CADD"
  },
  "variants": {
    "rs121913343": {
      "cadd_raw_score": 4.123,
      "cadd_phred_score": 25.67,
      "pathogenicity_level": "MEDIUM",
      "is_pathogenic": true,
      "chromosome": "7",
      "position": 140453136,
      "reference_allele": "A",
      "alternate_allele": "T",
      "gene": "BRAF",
      "protein_change": "p.Val600Glu"
    }
  }
}
```

## Configuration Options

### CADD Client Configuration

```python
config = {
    'use_cadd_api': True,           # Use CADD API (default: True)
    'cadd_cache_dir': 'cache/cadd', # Cache directory
    'batch_size': 100               # API batch size
}

scorer = PathogenicityScorer(config)
```

### Offline Mode

For environments without internet access:

```python
config = {
    'use_cadd_api': False,  # Use local CADD data only
    'cadd_cache_dir': 'cache/cadd'
}

scorer = PathogenicityScorer(config)
```

## Performance Optimization

### Caching

The CADD scorer automatically caches results in a SQLite database:

- **Location**: `cache/cadd/cadd_cache.db`
- **Benefits**: Faster subsequent queries for same variants
- **Automatic**: No manual intervention required

### Batch Processing

For large variant sets:

```python
# Process in batches (automatic)
results = scorer.score_variants(large_variant_list)

# Monitor progress
for batch_idx, batch in enumerate(variant_batches):
    print(f"Processing batch {batch_idx + 1}/{len(variant_batches)}")
```

## Error Handling

### Common Issues

1. **Network Connectivity**
   ```python
   # CADD API unavailable
   config = {'use_cadd_api': False}
   scorer = PathogenicityScorer(config)
   ```

2. **Invalid Variant Format**
   ```python
   # Ensure required fields
   variant = {
       'chromosome': '7',
       'position': 140453136,
       'reference_allele': 'A',
       'alternate_allele': 'T'
   }
   ```

3. **Missing CADD Data**
   ```python
   # Download pre-computed scores
   # See Installation section above
   ```

## Integration with Existing Pipeline

### With Variant Impact Analyzer

```python
from modules.integrated_pathogenicity_analyzer import IntegratedPathogenicityAnalyzer

analyzer = IntegratedPathogenicityAnalyzer()

# Combines CADD + existing analysis
results = analyzer.analyze_variants_comprehensive(variants)
```

### With Streamlit Viewer

```python
# Add CADD scores to existing results
import streamlit as st
from modules.cadd_pathogenicity_scorer import PathogenicityScorer

scorer = PathogenicityScorer()
cadd_results = scorer.score_variants(variants)

# Display in Streamlit
st.write("CADD Pathogenicity Scores")
for variant_id, score in cadd_results.items():
    st.metric(f"CADD Phred: {variant_id}", f"{score.cadd_phred_score:.2f}")
```

## Clinical Interpretation Guidelines

### High Pathogenicity (CADD Phred ≥ 30)
- **Action**: Prioritize for clinical review
- **Consider**: Functional studies, family segregation
- **Documentation**: Strong evidence for pathogenicity

### Medium Pathogenicity (CADD Phred 20-30)
- **Action**: Review in clinical context
- **Consider**: Additional evidence (segregation, functional studies)
- **Documentation**: Moderate evidence for pathogenicity

### Low Pathogenicity (CADD Phred 10-20)
- **Action**: Monitor for additional evidence
- **Consider**: Population frequency, conservation
- **Documentation**: Limited evidence

### Benign (CADD Phred < 10)
- **Action**: Likely benign, routine monitoring
- **Consider**: Population frequency
- **Documentation**: Evidence against pathogenicity

## Best Practices

1. **Always validate input variants** before scoring
2. **Use CADD API when possible** for most accurate results
3. **Cache results** for efficiency in batch processing
4. **Combine with other evidence** for comprehensive assessment
5. **Document confidence levels** in clinical reports
6. **Regular updates** of CADD data for latest annotations

## References

- **CADD Paper**: Kircher et al. (2014) Nature Genetics
- **CADD Website**: https://cadd.gs.washington.edu/
- **CADD API**: https://cadd.gs.washington.edu/api/v1.0/

## Support

For issues or questions:

1. Check the error logs in `logs/`
2. Verify network connectivity for CADD API
3. Ensure proper variant format
4. Contact: joshdrobert@gmail.com

---

**Note**: CADD is the gold standard for variant pathogenicity assessment and is widely used in clinical genomics. This implementation provides efficient, standardized access to CADD scoring with proper caching and error handling. 