import sqlite3
import re
from collections import Counter

DB_PATH = 'cache/alphamissense/alphamissense_hg38.db'

# Test BRCA1 var_002 variant
VARIANT_ID = "var_002"
GENE = "BRCA1"
PROTEIN_CHANGE = "p.Cys61Gly"
REFERENCE = "T"
ALTERNATE = "C"
VEP_TRANSCRIPT_IDS = [
    "ENST00000497488",
    "ENST00000489037",
    "ENST00000478531",
    "ENST00000357654",
    "ENST00000473961",
    "ENST00000477152",
    "ENST00000352993",
    "ENST00000493919",
    "ENST00000494123",
    "ENST00000471181",
    "ENST00000652672",
    "ENST00000634433",
    "ENST00000476777",
    "ENST00000700081",
    "ENST00000470026",
    "ENST00000713676",
    "ENST00000618469",
    "ENST00000461574",
    "ENST00000644555",
    "ENST00000468300",
    "ENST00000700082",
    "ENST00000644379",
    "ENST00000484087",
    "ENST00000586385",
    "ENST00000591534",
    "ENST00000591849",
    "ENST00000493795",
    "ENST00000461221",
    "ENST00000491747",
    "ENST00000472490",
    "ENST00000700182",
    "ENST00000621897",
    "ENST00000354071",
    "ENST00000700183",
    "ENST00000492859",
    "ENST00000642945",
    "ENST00000700184",
    "ENST00000461798",
    "ENST00000700083",
    "ENST00000700185",
    "ENST00000700186"
]

def test_brca1_var002_filtering():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # Extract position number from protein change
    position_match = re.search(r'(\d+)', PROTEIN_CHANGE)
    if not position_match:
        print(f"Could not extract position from protein change: {PROTEIN_CHANGE}")
        return
    
    position_number = position_match.group(1)
    position_pattern = f"%{position_number}%"
    
    print(f"Testing BRCA1 variant: {PROTEIN_CHANGE}")
    print(f"Position number: {position_number}")
    print(f"Position pattern: {position_pattern}")
    print(f"Reference: {REFERENCE} | Alternate: {ALTERNATE}")
    print("-" * 80)
    
    # Create placeholders for the IN clause
    placeholders = ','.join(['?' for _ in VEP_TRANSCRIPT_IDS])
    
    # Initial query - position-based matching
    query = f"""
        SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
        FROM alphamissense 
        WHERE transcript_id IN ({placeholders})
        AND protein_variant LIKE ?
        ORDER BY am_pathogenicity DESC
    """
    
    params = VEP_TRANSCRIPT_IDS + [position_pattern]
    cursor.execute(query, params)
    matches = cursor.fetchall()
    
    print(f"Initial position-based matches: {len(matches)}")
    
    # If too many matches, apply filtering
    if matches and len(matches) > 5:
        print(f"Too many matches ({len(matches)}), applying filtering...")
        
        # Try strict protein variant matching
        strict_query = f"""
            SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
            FROM alphamissense 
            WHERE transcript_id IN ({placeholders})
            AND protein_variant = ?
            ORDER BY am_pathogenicity DESC
        """
        strict_params = VEP_TRANSCRIPT_IDS + [PROTEIN_CHANGE]
        cursor.execute(strict_query, strict_params)
        strict_matches = cursor.fetchall()
        
        if strict_matches:
            matches = strict_matches
            print(f"Strict matching successful: {len(matches)} matches")
        else:
            print("Strict matching failed, trying REF/ALT filtering...")
            # Try filtering by REF and ALT
            ref_alt_query = f"""
                SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                FROM alphamissense 
                WHERE transcript_id IN ({placeholders})
                AND protein_variant LIKE ?
                AND REF = ?
                AND ALT = ?
                ORDER BY am_pathogenicity DESC
            """
            ref_alt_params = VEP_TRANSCRIPT_IDS + [position_pattern, REFERENCE, ALTERNATE]
            cursor.execute(ref_alt_query, ref_alt_params)
            ref_alt_matches = cursor.fetchall()
            if ref_alt_matches:
                if len(ref_alt_matches) > 10:
                    print(f"REF/ALT filtering successful but too many matches ({len(ref_alt_matches)}), trying pathogenicity filtering...")
                    # Filter by pathogenicity score
                    pathogenicity_query = f"""
                        SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                        FROM alphamissense 
                        WHERE transcript_id IN ({placeholders})
                        AND protein_variant LIKE ?
                        AND REF = ?
                        AND ALT = ?
                        AND am_pathogenicity > 0.8
                        ORDER BY am_pathogenicity DESC
                        LIMIT 10
                    """
                    pathogenicity_params = VEP_TRANSCRIPT_IDS + [position_pattern, REFERENCE, ALTERNATE]
                    cursor.execute(pathogenicity_query, pathogenicity_params)
                    pathogenicity_matches = cursor.fetchall()
                    if pathogenicity_matches:
                        matches = pathogenicity_matches
                        print(f"Pathogenicity filtering successful: {len(matches)} matches (score > 0.8)")
                    else:
                        matches = ref_alt_matches
                        print(f"Pathogenicity filtering failed, using REF/ALT matches: {len(matches)} matches")
                else:
                    matches = ref_alt_matches
                    print(f"REF/ALT filtering successful: {len(matches)} matches")
            else:
                print("REF/ALT filtering failed, trying pathogenicity filtering...")
                # Filter by pathogenicity score
                pathogenicity_query = f"""
                    SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant, REF, ALT
                    FROM alphamissense 
                    WHERE transcript_id IN ({placeholders})
                    AND protein_variant LIKE ?
                    AND am_pathogenicity > 0.8
                    ORDER BY am_pathogenicity DESC
                    LIMIT 10
                """
                pathogenicity_params = VEP_TRANSCRIPT_IDS + [position_pattern]
                cursor.execute(pathogenicity_query, pathogenicity_params)
                pathogenicity_matches = cursor.fetchall()
                if pathogenicity_matches:
                    matches = pathogenicity_matches
                    print(f"Pathogenicity filtering successful: {len(matches)} matches (score > 0.8)")
                else:
                    print("Pathogenicity filtering failed, using original broad matches")
    else:
        print(f"Matches within acceptable range ({len(matches)}), no filtering needed")
    
    # Display results
    print(f"\nFinal results: {len(matches)} matches")
    print("-" * 80)
    
    if matches:
        # Calculate statistics
        pathogenicity_scores = [match[0] for match in matches if match[0] is not None]
        classes = [match[1] for match in matches if match[1] is not None]
        matching_transcripts = list(set([match[3] for match in matches]))
        
        if pathogenicity_scores:
            avg_score = sum(pathogenicity_scores) / len(pathogenicity_scores)
            print(f"Average pathogenicity score: {avg_score:.5f}")
        
        if classes:
            class_counter = Counter(classes)
            most_common_class = class_counter.most_common(1)[0][0]
            print(f"Most common class: {most_common_class}")
        
        print(f"Matching transcripts: {matching_transcripts}")
        
        # Show all results
        print("\nAll results:")
        for i, match in enumerate(matches, 1):
            score, class_val, uniprot, transcript, protein_var, ref, alt = match
            print(f"  {i}. {protein_var} | Score: {score:.4f} | Class: {class_val} | Transcript: {transcript} | REF: {ref} | ALT: {alt}")
    else:
        print("No matches found")
    
    conn.close()

if __name__ == "__main__":
    test_brca1_var002_filtering() 