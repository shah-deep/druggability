import streamlit as st
import json

# Load results
with open("full_results.json") as f:
    results = json.load(f)

st.title("🧪 Mechanistic Coherence Viewer")

# Summary
st.header("Final Analysis")
st.metric("Final Score", results["combined_score"]["combined_score"])
st.metric("Mechanistic Score", results["combined_score"]["mechanistic_score"])
st.metric("Best Druggability Score", results["combined_score"]["best_druggability_score"])
st.write(f"**Recommendation:** {results['recommendation']}")

# Mechanism Summary
summary = results["compound_analysis"]["mechanistic_coherence"].get("summary", "")
st.header("🧬 Mechanism Explanation")
if summary:
    for line in summary.split('. '):
        line = line.strip()
        if line:
            st.markdown(f"- {line.rstrip('.')}.")
else:
    st.info("No mechanism summary available.")

# Category scores
st.header("Category Scores")
for cat, score in results["compound_analysis"]["mechanistic_coherence"]["category_scores"].items():
    st.write(f"- **{cat.title()}**: {score:.3f}")

# Conflicts
st.header("Conflict Flags")
conflicts = results["compound_analysis"]["mechanistic_coherence"].get("conflict_flags", [])
if conflicts:
    for c in conflicts:
        st.warning(f"⚠️ {c}")
else:
    st.success("✅ No conflicts detected")

# Model Outputs
st.header("Model Outputs")
for cat, outputs in results["compound_analysis"]["mechanistic_coherence"]["model_outputs"].items():
    st.subheader(f"{cat.title()} Models")
    st.json(outputs)

# Druggability
if "druggability_assessment" in results["compound_analysis"]:
    st.header("Druggability Analysis (Top Pockets)")
    for pid, pocket in list(results["compound_analysis"]["druggability_assessment"].items())[:3]:
        st.subheader(f"Pocket {pid}")
        st.json(pocket)
