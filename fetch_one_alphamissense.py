import sqlite3
import re

DB_PATH = 'cache/alphamissense/alphamissense_hg38.db'

# Test CFTR variant specifically
VARIANT_ID = "var_CFTR_DF508"
GENE = "CFTR"
PROTEIN_CHANGE = "p.Phe508del"
VEP_TRANSCRIPT_IDS = [
    "ENST00000673785", "ENST00000436097", "ENST00000546407", "ENST00000446805", "ENST00000699596",
    "ENST00000699597", "ENST00000699598", "ENST00000699599", "ENST00000647720", "ENST00000699600",
    "ENST00000699585", "ENST00000699601", "ENST00000699602", "ENST00000685018", "ENST00000687278",
    "ENST00000693465", "ENST00000699603", "ENST00000693480", "ENST00000692802", "ENST00000647639",
    "ENST00000649850", "ENST00000648260", "ENST00000649406", "ENST00000647978", "ENST00000649781",
    "ENST00000699604", "ENST00000699605", "ENST00000003084", "ENST00000426809", "ENST00000472848",
    "ENST00000468795", "ENST00000689011", "ENST00000699606", "ENST00000600166", "ENST00000429014",
    "ENST00000608965", "ENST00000610149", "ENST00000621535"
]

def test_cftr_filtering():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # Extract position number from protein change
    position_match = re.search(r'(\d+)', PROTEIN_CHANGE)
    if not position_match:
        print(f"Could not extract position from protein change: {PROTEIN_CHANGE}")
        return
    
    position_number = position_match.group(1)
    position_pattern = f"%{position_number}%"
    
    print(f"Testing CFTR variant: {PROTEIN_CHANGE}")
    print(f"Position number: {position_number}")
    print(f"Position pattern: {position_pattern}")
    print("-" * 80)
    
    # Create placeholders for the IN clause
    placeholders = ','.join(['?' for _ in VEP_TRANSCRIPT_IDS])
    
    # Initial query - position-based matching
    query = f"""
        SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant
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
            SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant
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
            print("Strict matching failed, trying pathogenicity filtering...")
            
            # Filter by pathogenicity score
            pathogenicity_query = f"""
                SELECT am_pathogenicity, am_class, uniprot_id, transcript_id, protein_variant
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
            from collections import Counter
            class_counter = Counter(classes)
            most_common_class = class_counter.most_common(1)[0][0]
            print(f"Most common class: {most_common_class}")
        
        print(f"Matching transcripts: {matching_transcripts}")
        
        # Show all results
        print("\nAll results:")
        for i, match in enumerate(matches, 1):
            score, class_val, uniprot, transcript, protein_var = match
            print(f"  {i}. {protein_var} | Score: {score:.4f} | Class: {class_val} | Transcript: {transcript}")
    else:
        print("No matches found")
    
    conn.close()

if __name__ == "__main__":
    test_cftr_filtering() 