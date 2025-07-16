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
   git clone https://github.com/shah-deep/druggability
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
   fpocket -h  # Verify fpocket is installed
   ```

4. **Set up environment variables (optional but recommended):**
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