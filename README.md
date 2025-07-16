# Enhanced Mechanistic Coherence Module

A comprehensive bioinformatics pipeline for analyzing genetic variants, protein structures, and pathway dynamics to assess drug target druggability and biological coherence across multiple scales.

## Project Overview

This project implements an enhanced mechanistic coherence analysis system that integrates:
- **Variant Impact Analysis**: Assessment of genetic variant effects on protein function
- **Structural Analysis**: Protein structure analysis and pocket detection
- **Pathway Dynamics**: Analysis of pathway-level impacts and network effects
- **Coherence Scoring**: Cross-scale biological consistency validation
- **Clinical Integration**: Tissue-specific biomarker relevance assessment

## Project Structure

```
druggability/
├── modules/                          # Core analysis modules
│   ├── enhanced_input_processor.py   # Input validation and preprocessing
│   ├── variant_impact_analyzer.py    # Variant effect prediction
│   ├── structure_function_integrator.py # Structural analysis integration
│   ├── sequence_variant_applier.py   # Sequence variant application
│   ├── coherence_analyzer.py         # Cross-scale coherence analysis (includes clinical_translator functionality)
│   ├── pathway_dynamics_analyzer.py  # Pathway-level impact analysis
│   ├── intelligent_coherence_scorer.py # Final scoring and integration
│   ├── pocket_detection.py           # Protein pocket identification
│   ├── feature_extraction.py         # Feature extraction utilities
│   ├── clinvar_annotator.py         # ClinVar database integration
│   ├── transcript_resolver.py        # Transcript mapping utilities
│   └── logging_config.py            # Logging configuration
├── orchestrator/                     # Pipeline orchestration
│   └── enhanced_pipeline_orchestrator.py # Main pipeline coordinator
├── examples/                         # Example input files
│   ├── ex_variants.json             # Example genetic variants
│   ├── ex_clinical_data.json        # Example clinical data
│   ├── ex_protein_seq.txt           # Example protein sequence
│   └── protein_42.pdb               # Example PDB structure
├── outputs/                          # Pipeline output files
├── tests/                           # Test suite
│   ├── test_variant_impact_analyzer.py
│   ├── test_structure_function_integrator.py
│   └── test_pathway_validation.py
├── ProteinMPNN/                     # ProteinMPNN integration
├── mech-coherence-env.yml           # Conda environment specification
└── README.md                        # This file
```

## Setup and Installation

### Prerequisites

- **Conda/Miniconda** - For environment management
- **Python 3.10** - Required Python version
- **Git** - For cloning the repository

### Installation Steps

1. **Clone the repository and navigate to the project directory:**
   ```bash
   git clone <repository-url>
   cd druggability
   ```

2. **Create and activate the conda environment:**
   ```bash
   conda env create -f mech-coherence-env.yml
   conda activate mech-coherence
   ```

3. **Verify the installation:**
   ```bash
   python -c "import numpy, pandas, torch, biopython; print('All core dependencies installed successfully')"
   ```

4. **Install fpocket (required for structural analysis):**
   ```bash
   conda install -c conda-forge fpocket -y
   ```

5. **Set up environment variables (optional but recommended):**
   Create a `.env` file in the project root:
   ```bash
   # Optional but recommended for enhanced functionality
   NCBI_API_KEY=your_ncbi_api_key_here
   ```

## Usage

### Running the Complete Pipeline

The main entry point is the `EnhancedPipelineOrchestrator`. Run the pipeline using:

```bash
# Activate the conda environment first
conda activate mech-coherence

# Run the pipeline
python -m orchestrator.enhanced_pipeline_orchestrator \
    --variants examples/ex_variants.json \
    --clinical examples/ex_clinical_data.json \
    --protein-seq examples/ex_protein_seq.txt \
    --pdb examples/protein_42.pdb \
    --output-dir outputs \
    --max-workers 8
```

### Input File Formats

#### Variants JSON (`ex_variants.json`)
```json
[
  {
    "id": "var_001",
    "position": 17498,
    "reference": "G",
    "alternate": "A",
    "gene": "TP53",
    "protein_change": "p.Arg175His"
  }
]
```

#### Clinical Data JSON (`ex_clinical_data.json`)
```json
{
  "patient_id": "P001",
  "age": 55,
  "gender": "female",
  "diagnosis": "breast_cancer",
  "stage": "II",
  "biomarkers": {
    "HER2": "positive",
    "ER": "positive",
    "PR": "negative"
  },
  "prior_treatments": ["chemotherapy", "radiation"]
}
```

#### Protein Sequence (`ex_protein_seq.txt`)
Plain text file containing the protein sequence in FASTA format.

#### PDB Structure (`protein_42.pdb`)
Standard PDB format file containing the protein structure.

### Output Files

The pipeline generates several output files in the specified output directory:
- `processed_input_YYYYMMDD_HHMMSS.json` - Processed and validated input data
- `variant_impact_results_YYYYMMDD_HHMMSS.json` - Variant impact analysis results
- `structural_results_YYYYMMDD_HHMMSS.json` - Structural analysis results
- `sequence_variant_results_YYYYMMDD_HHMMSS.json` - Sequence variant analysis
- `coherence_results_YYYYMMDD_HHMMSS.json` - Cross-scale coherence analysis
- `pathway_dynamics_results_YYYYMMDD_HHMMSS.json` - Pathway dynamics analysis
- `final_results_YYYYMMDD_HHMMSS.json` - Final integrated results

## Important Notes

### Large File Dependencies

⚠️ **AlphaMissense Database**: The current setup requires the AlphaMissense database file:
```
cache/alphamissense/alphamissense_hg38.db
```

This file is **12GB** and too large to distribute with the repository. For demonstration purposes, please contact the Deep Shah for a live demo via video call.

### Module Integration

- **Clinical Translator**: The functionality of `clinical_translator.py` has been integrated into `coherence_analyzer.py` to streamline the analysis pipeline.

### Testing

Run the test suite using pytest:

```bash
# Activate the conda environment first
conda activate mech-coherence

# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_variant_impact_analyzer.py

# Run with verbose output
pytest -v tests/
```

The test suite can be integrated into CI/CD pipelines for automated testing.

## Troubleshooting

### Common Issues

1. **Conda environment creation fails**: 
   - Make sure you have the latest version of conda installed
   - If you see "Solving environment" taking too long, the environment file has been optimized to avoid dependency conflicts

2. **Import errors**: 
   - Ensure you're using the correct conda environment (`conda activate mech-coherence`)
   - Verify installation with: `python -c "import numpy, pandas, torch, biopython; print('OK')"`

3. **Memory issues**: 
   - Reduce `--max-workers` parameter for systems with limited RAM

4. **Missing fpocket**: 
   - Install fpocket using: `conda install -c conda-forge fpocket -y`
   - Verify installation with: `fpocket -h`

5. **gseapy installation issues**: 
   - The environment file excludes `gseapy` due to Rust compiler requirements
   - If you need pathway enrichment analysis, install Rust first: `brew install rust` (macOS) or `conda install rust`
   - Then install gseapy: `pip install gseapy`

6. **fpocket error**: 
   - If you see "No such file or directory: 'fpocket'", install fpocket using the command above
   - Make sure you're in the conda environment when running the pipeline

### System Requirements

- **Minimum RAM**: 8GB (16GB recommended for large datasets)
- **Storage**: At least 20GB free space for outputs and temporary files
- **CPU**: Multi-core processor recommended for parallel processing
- **Operating System**: macOS, Linux, or Windows with WSL

### Environment Details

The conda environment includes:
- **Core scientific packages**: numpy, pandas, scipy, scikit-learn
- **Bioinformatics tools**: biopython, scanpy, anndata
- **Machine learning**: torch, transformers
- **Visualization**: matplotlib, seaborn
- **Web frameworks**: streamlit, aiohttp
- **Structural biology**: mdtraj, cobra, fpocket
- **Testing**: pytest

All packages are configured for compatibility and should install without conflicts.