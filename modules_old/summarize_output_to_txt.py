import json
import requests
from pathlib import Path

# Config
API_URL = "http://localhost:1234/v1/chat/completions"  # LM Studio default
MODEL_NAME = "local-model"  # You can leave it or change it to your exact LM Studio model

# File paths
project_dir = Path(__file__).parent
pockets_file = project_dir / "pockets_protein_42.json"
features_file = project_dir / "pocket_features_protein_42.json"
scores_file = project_dir / "druggability_scores_protein_42.json"
output_txt = project_dir / "druggability_summary.txt"

# Prompt Template
prompt = f"""
You are an expert in structural bioinformatics and drug discovery.

You are reviewing the results of a binding pocket druggability pipeline applied to the protein `protein_42.pdb`.

You have access to the following:
1. `pockets_protein_42.json` – pocket locations and sizes
2. `pocket_features_protein_42.json` – features like volume, shape complexity, hydrophobicity, polarity, and solvent accessibility
3. `druggability_scores_protein_42.json` – final druggability scores for each pocket

Your task is to write a **detailed, professional, and assertive analysis**. Avoid vague commentary or speculation. Focus on clear findings, rankings, and practical implications for drug design. Be concise, confident, and specific.

**Include:**
- Total number of pockets found.
- Top 3 pockets with the highest scores (include scores).
- How pocket features support or undermine their druggability.
- Which pocket(s) should be prioritized for structure-based drug design and why.
- Any red flags or pockets that should be avoided due to poor properties.
- Final recommendation for next steps (e.g., virtual screening, modeling, etc.).

Write for a technically literate but non-expert audience (e.g., a biotech project manager or interdisciplinary research team).
"""

# Build messages
messages = [
    {"role": "system", "content": "You are a helpful assistant."},
    {"role": "user", "content": prompt}
]

# Call LM Studio (OpenAI-compatible API)
response = requests.post(API_URL, json={
    "model": MODEL_NAME,
    "messages": messages,
    "temperature": 0.7
})

if response.status_code == 200:
    reply = response.json()["choices"][0]["message"]["content"]
    output_txt.write_text(reply)
    print(f"✅ Summary saved to {output_txt}")
else:
    print(f"❌ Error from LM Studio: {response.status_code}")
    print(response.text)