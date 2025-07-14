# ClinVar Annotation Module

This module provides comprehensive ClinVar variant annotation functionality for your variant analysis pipeline.

## Overview

The ClinVar annotation module supports multiple lookup strategies to find clinical significance information for genetic variants:

1. **Known Variants Cache** - Pre-defined pathogenic variants
2. **Local Database Cache** - SQLite cache of previously queried variants
3. **Local Summary File** - ClinVar variant summary file (if available)
4. **ClinVar API** - NCBI E-utilities API for real-time queries
5. **Fallback Mechanisms** - Graceful handling when no annotation is found

## Features

- **Multiple Lookup Strategies**: Tries different sources in order of reliability
- **Caching**: Local SQLite database for performance
- **API Integration**: NCBI E-utilities for real-time ClinVar queries
- **Comprehensive Results**: Clinical significance, review status, conditions, etc.
- **Error Handling**: Graceful fallbacks and detailed warnings
- **Batch Processing**: Annotate multiple variants efficiently

## Installation

The module requires the following Python packages:

```bash
pip install requests
```

## Usage

### Basic Usage

```python
from modules.clinvar_annotator import ClinVarAnnotator

# Initialize annotator
annotator = ClinVarAnnotator()

# Your variant data
variant = {
    "id": "var_001",
    "position": 7675088,
    "reference": "G",
    "alternate": "A",
    "gene": "TP53",
    "protein_change": "p.Arg175His",
    "chromosome": "17",
    "transcript_id": "ENST00000269305"
}

# Annotate single variant
result = annotator.annotate_variant(variant)

print(f"Clinical Significance: {result.clinical_significance}")
print(f"Review Status: {result.review_status}")
print(f"Condition: {result.condition}")
print(f"Source: {result.source}")
```

### Batch Processing

```python
# Annotate multiple variants
variants = [
    {
        "id": "var_001",
        "gene": "TP53",
        "protein_change": "p.Arg175His",
        # ... other fields
    },
    {
        "id": "var_002", 
        "gene": "CFTR",
        "protein_change": "p.Phe508del",
        # ... other fields
    }
]

results = annotator.annotate_variants(variants)

for variant_id, result in results.items():
    print(f"{variant_id}: {result.clinical_significance}")
```

## Input Format

The annotator expects variant dictionaries with the following fields:

### Required Fields

- `gene`: Gene symbol (e.g., "TP53")
- `protein_change`: Protein change in HGVS format (e.g., "p.Arg175His")

### Optional Fields

- `id`: Variant identifier
- `position`: Genomic position
- `reference`: Reference allele
- `alternate`: Alternate allele
- `chromosome`: Chromosome number
- `transcript_id`: Transcript identifier

## Output Format

The annotator returns `ClinVarAnnotation` objects with the following fields:

### Core Fields

- `variant_id`: Variant identifier
- `gene`: Gene symbol
- `protein_change`: Protein change
- `clinical_significance`: Clinical significance (Pathogenic, Likely_pathogenic, etc.)
- `review_status`: Review status from ClinVar
- `condition`: Associated condition/disease
- `variation_id`: ClinVar variation ID
- `source`: Data source (known_variants, local_cache, api, etc.)

### Additional Fields

- `last_evaluated`: Date last evaluated
- `submitter`: Submitter information
- `star_rating`: Review star rating
- `accession`: ClinVar accession number
- `warnings`: List of warnings/errors

## Example Results

For the TP53 p.Arg175His variant:

```python
# Example output
ClinVarAnnotation(
    variant_id="var_001",
    gene="TP53",
    protein_change="p.Arg175His",
    clinical_significance="Pathogenic",
    review_status="criteria provided, multiple submitters, no conflicts",
    condition="Li-Fraumeni syndrome 1",
    variation_id="376649",
    star_rating=4,
    source="known_variants"
)
```

## Configuration

### Environment Variables

Set your NCBI API key for better rate limits:

```bash
export NCBI_API_KEY="your_api_key_here"
```

### Cache Configuration

The module uses a local SQLite database for caching. Cache files are stored in:

- Database: `cache/clinvar/clinvar_cache.db`
- Summary file: `cache/clinvar_variant_summary.txt.gz` (optional)

### Custom Configuration

```python
# Custom cache directory and summary file
annotator = ClinVarAnnotator(
    cache_dir="my_cache/clinvar",
    summary_file="my_data/clinvar_summary.txt.gz"
)
```

## API Integration

The module integrates with NCBI E-utilities for real-time ClinVar queries:

- **Search**: Uses `esearch.fcgi` to find variant IDs
- **Fetch**: Uses `efetch.fcgi` to get detailed variant information
- **Rate Limiting**: Built-in delays to respect API limits
- **Error Handling**: Graceful fallbacks when API is unavailable

## Performance Considerations

1. **Caching**: Results are cached in SQLite database for fast subsequent lookups
2. **Batch Processing**: Use `annotate_variants()` for multiple variants
3. **API Limits**: Respect NCBI rate limits (set API key for higher limits)
4. **Local Data**: Use local ClinVar summary file for offline operation

## Error Handling

The module provides comprehensive error handling:

- **Network Errors**: Graceful fallback to cached/local data
- **API Errors**: Detailed logging and fallback strategies
- **Missing Data**: Clear warnings when no annotation is found
- **Invalid Input**: Validation and helpful error messages

## Integration with Existing Pipeline

The ClinVar annotator can be easily integrated into your existing variant analysis pipeline:

```python
# In your main analysis pipeline
from modules.clinvar_annotator import ClinVarAnnotator

def analyze_variants(variants):
    # Initialize annotators
    clinvar_annotator = ClinVarAnnotator()
  
    # Get ClinVar annotations
    clinvar_results = clinvar_annotator.annotate_variants(variants)
  
    # Combine with other analysis results
    for variant in variants:
        variant_id = variant['id']
        clinvar_annotation = clinvar_results[variant_id]
      
        # Add ClinVar data to your results
        variant['clinvar_significance'] = clinvar_annotation.clinical_significance
        variant['clinvar_review_status'] = clinvar_annotation.review_status
        variant['clinvar_condition'] = clinvar_annotation.condition
  
    return variants
```

## Testing

Run the test scripts to verify functionality:

```bash
# Test basic functionality
python test_clinvar_annotation.py

# Test integration example
python example_clinvar_usage.py
```

## Troubleshooting

### Common Issues

1. **No annotations found**: Check if gene and protein_change are correctly formatted
2. **API errors**: Verify network connection and API key (if using)
3. **Cache issues**: Delete cache database to reset: `rm cache/clinvar/clinvar_cache.db`
4. **Rate limiting**: Set NCBI API key for higher limits

### Debug Mode

Enable detailed logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Dependencies

- `requests`: HTTP requests for API calls
- `sqlite3`: Local caching (built-in)
- `xml.etree.ElementTree`: XML parsing (built-in)
- `gzip`: Compressed file handling (built-in)

## License

This module is part of the druggability analysis pipeline and follows the same license terms.
