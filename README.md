# Mechanistic Coherence & Druggability Scoring Pipeline

This repository contains a modular Python pipeline for evaluating drug candidates through two key components:

- **Druggability Assessment**: Detects and scores potential binding pockets using structural features extracted from protein `.pdb` files.
- **Mechanistic Coherence Module**: Integrates multi-omic and structural models to assess biological plausibility across genotype, structure, simulation, biomarker, and metabolic dimensions.

The pipeline outputs structured JSON scores, conflict flags, summaries, and recommendations to support early-stage therapeutic evaluation.

---

## 🚫 Proprietary License — Round 2 Evaluation Only

Copyright (c) 2025 Joshua Robert  
All rights reserved.

This source code and all associated files are the intellectual property of Joshua Robert.  
It is provided **solely for evaluation and review purposes** by Convexia.  
Use of this code for **any commercial, research, or production purposes** is strictly prohibited without the express written consent of the author.

By accessing this repository, you agree:
- Not to copy, redistribute, or modify this code beyond the scope of evaluation.
- Not to incorporate this code into any proprietary system or product without permission.
- That all IP rights remain with the author unless otherwise agreed in writing.

To obtain a commercial, research, or collaborative license, please contact:  
📧 joshdrobert@gmail.com

---

## ⚙️ 1. Setup Instructions

### ✅ Python Version

- Python **3.10** is required.

### ✅ Install Dependencies

```bash
# Clone the repository and enter the directory
git clone <your-repo-url>
cd <your-repo-folder>

# Create and activate a virtual environment
python3 -m venv venv
source venv/bin/activate        # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

-- 

## External Tools (Required Separately)
Ensure the following are installed and accessible from your system PATH:

RDKit – for molecular descriptors

OpenBabel – for file conversion (optional but helpful)

PyMOL – for visualizing protein-ligand complexes

🚀 2. Execution Instructions
Run the pipeline step-by-step using your .pdb file (e.g., protein_42.pdb):

1. Detect Binding Pockets
python pocket_detection.py protein_42.pdb -o protein_42_out/

2. Extract Features for Druggability Scoring
python feature_extraction.py protein_42_out/pockets_protein_42.json protein_42.pdb -d protein_42_out/pockets/

3. Score Pockets
python scoring_engine.py protein_42_out/pocket_features_protein_42.json

4. Evaluate Mechanistic Coherence
python its42_mechanistic_coherence.py protein_42.pdb \
    --dna-sequence "ATGCGTACGTACGTACGTACGTACG..." \
    --rna-sequence "AUGCGUACGUACGUACGUACGUACG..." \
    --protein-sequence "MKTAYIAKQRQISFVKSHFSRQLE..." \
    --clinical-data ex_clinical_data.json \
    --missense-variants ex_variants.json \
    -o full_results.json

5. Launch Results Viewer
streamlit run results_viewer.py

📂 3. Inputs
File Type	Description	Required?
*.pdb	Protein structure file	✅ Required
*.json	Pocket detection output (pockets_*.json)	✅ Required
*.json	Clinical data and variant files for ITS4.2 module	Optional
Sequences	DNA / RNA / protein sequences (as strings)	✅ Required

📁 4. Outputs
Output File	Description	Location
pockets_protein_42.json	Binding pocket detection output	protein_42_out/
pocket_features_protein_42.json	Feature vectors for each pocket	protein_42_out/
scored_pockets.json	Druggability scores with ranking	protein_42_out/
full_results.json	Full mechanistic coherence report (ITS4.2)	User-specified
Streamlit UI	Interactive results browser	Opens in browser

🛠️ 5. Troubleshooting
❌ ModuleNotFoundError: No module named 'rdkit'
RDKit is not available via pip. You must install it separately (e.g., via conda):

bash
Copy
Edit
conda install -c rdkit rdkit
❌ streamlit: command not found
You’re likely outside your virtual environment. Activate it and install:

bash
Copy
Edit
source venv/bin/activate      # or venv\Scripts\activate
pip install streamlit
❌ FileNotFoundError: pockets_protein_42.json not found
Make sure you run pocket_detection.py first, and verify your paths.

❌ Invalid or empty sequence
Check that your --dna-sequence, --rna-sequence, and --protein-sequence values are properly formatted and non-empty strings.

If you encounter any other issues, feel free to contact the author for clarification.

-- 

📫 Contact
Joshua Robert
📧 joshdrobert@gmail.com
🔬 M.D. / M.Eng. candidate | Texas A&M EnMed
💡 AI + Biotech Innovator
