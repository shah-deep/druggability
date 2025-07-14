import pandas as pd
import sqlite3
import os

def check_database_status():
    """
    Check the current status of the AlphaMissense database
    """
    db_file = "cache/alphamissense/alphamissense_hg38.db"
    
    if not os.path.exists(db_file):
        print("Database file does not exist yet.")
        return None
    print("Connecting to database...")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    print("Checking if table exists...")
    # Check if table exists
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='alphamissense'")
    if not cursor.fetchone():
        print("Database exists but table 'alphamissense' not found.")
        conn.close()
        return None
    print("Getting current row count...")
    # Get current row count
    cursor.execute("SELECT COUNT(*) FROM alphamissense")
    current_rows = cursor.fetchone()[0]
    
    # Check if indexes exist
    cursor.execute("SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'")
    indexes = [row[0] for row in cursor.fetchall()]
    
    print(f"Current database status:")
    print(f"  - Total rows: {current_rows:,}")
    print(f"  - Indexes: {indexes}")
    
    # Create indexes if they don't exist
    if not indexes:
        print("\nCreating indexes...")
        indexes_to_create = [
            ("idx_transcript_id", "transcript_id"),
            ("idx_protein_variant", "protein_variant"),
            ("idx_uniprot_id", "uniprot_id")
        ]
        
        for index_name, column in indexes_to_create:
            try:
                conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON alphamissense({column})")
                print(f"Created index: {index_name}")
            except Exception as e:
                print(f"Error creating index {index_name}: {e}")
        
        conn.commit()
    
    conn.close()
    return db_file

def search_alphamissense(transcript_ids, protein_variant_pattern=None, limit=100):
    """
    Search AlphaMissense database for specific transcripts and optionally protein variants
    """
    db_file = "cache/alphamissense/alphamissense_hg38.db"
    
    if not os.path.exists(db_file):
        print(f"Database file not found: {db_file}")
        return None
    
    conn = sqlite3.connect(db_file)
    
    # Build query
    if protein_variant_pattern:
        query = """
        SELECT * FROM alphamissense 
        WHERE transcript_id IN ({}) 
        AND protein_variant LIKE ?
        LIMIT ?
        """.format(','.join(['?' for _ in transcript_ids]))
        params = transcript_ids + [protein_variant_pattern, limit]
    else:
        query = """
        SELECT * FROM alphamissense 
        WHERE transcript_id IN ({}) 
        LIMIT ?
        """.format(','.join(['?' for _ in transcript_ids]))
        params = transcript_ids + [limit]
    
    result = pd.read_sql_query(query, conn, params=params)
    conn.close()
    
    return result

def test_search_function():
    """
    Test the search function with CFTR transcripts and F508
    """
    print("\n" + "="*50)
    print("TESTING SEARCH FUNCTION")
    print("="*50)
    
    cftr_transcripts = [
        "ENST00000003084", "ENST00000426809", "ENST00000441019", "ENST00000472848",
        "ENST00000647720", "ENST00000647978", "ENST00000648260", "ENST00000649406",
        "ENST00000649781", "ENST00000685018", "ENST00000687278", "ENST00000699585",
        "ENST00000699596", "ENST00000699597", "ENST00000699598", "ENST00000699599",
        "ENST00000699600", "ENST00000699601", "ENST00000699602", "ENST00000699604",
        "ENST00000699605"
    ]
    
    # Search for F508 variants
    print("Searching for F508 variants in CFTR transcripts...")
    results = search_alphamissense(cftr_transcripts, protein_variant_pattern="%F508%")
    if results is not None and len(results) > 0:
        print(f"Found {len(results)} F508 variants:")
        print(results)
    else:
        print("No F508 variants found in CFTR transcripts")
    
    # Also search for any variants in CFTR transcripts
    print("\nSearching for any variants in CFTR transcripts...")
    results = search_alphamissense(cftr_transcripts, limit=10)
    if results is not None and len(results) > 0:
        print(f"Found {len(results)} variants in CFTR transcripts:")
        print(results)
    else:
        print("No variants found in CFTR transcripts")

if __name__ == "__main__":
    # Check current status
    print("Checking current database status...")
    db_file = check_database_status()
    
    if db_file:
        # Test the search function
        test_search_function()
    else:
        print("Database not ready for searching.") 