import pandas as pd
import sqlite3
import os
from pathlib import Path

def create_alphamissense_database():
    """
    Convert AlphaMissense_hg38.tsv to SQLite database with indexes for fast searching
    """
    
    # File paths
    tsv_file = "cache/alphamissense/AlphaMissense_hg38.tsv"
    db_file = "cache/alphamissense/alphamissense_hg38.db"
    
    print(f"Reading TSV file: {tsv_file}")
    
    # Read the TSV file with proper column names
    # Skip the copyright header lines and use the actual column names
    df = pd.read_csv(tsv_file, 
                     sep='\t', 
                     comment='#',
                     names=['CHROM', 'POS', 'REF', 'ALT', 'genome', 'uniprot_id', 
                           'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class'])
    
    print(f"Loaded {len(df)} rows")
    print(f"Columns: {list(df.columns)}")
    print(f"Sample data:")
    print(df.head())
    
    # Create database connection
    conn = sqlite3.connect(db_file)
    
    # Write to SQLite database
    df.to_sql('alphamissense', conn, if_exists='replace', index=False)
    
    # Create indexes for fast searching
    print("Creating indexes...")
    
    # Index on transcript_id for fast transcript searches
    conn.execute("CREATE INDEX IF NOT EXISTS idx_transcript_id ON alphamissense(transcript_id)")
    
    # Index on protein_variant for fast protein variant searches
    conn.execute("CREATE INDEX IF NOT EXISTS idx_protein_variant ON alphamissense(protein_variant)")
    
    # Index on uniprot_id for protein searches
    conn.execute("CREATE INDEX IF NOT EXISTS idx_uniprot_id ON alphamissense(uniprot_id)")
    
    # Note: idx_chrom_pos and idx_am_class indexes removed as requested
    
    # Commit changes
    conn.commit()
    
    # Test the database
    print("\nTesting database...")
    
    # Test transcript search (similar to your PowerShell command)
    test_transcripts = [
        "ENST00000003084", "ENST00000426809", "ENST00000441019", "ENST00000472848",
        "ENST00000647720", "ENST00000647978", "ENST00000648260", "ENST00000649406",
        "ENST00000649781", "ENST00000685018", "ENST00000687278", "ENST00000699585",
        "ENST00000699596", "ENST00000699597", "ENST00000699598", "ENST00000699599",
        "ENST00000699600", "ENST00000699601", "ENST00000699602", "ENST00000699604",
        "ENST00000699605"
    ]
    
    # Search for F508 (Phe508del) in CFTR transcripts
    query = """
    SELECT * FROM alphamissense 
    WHERE transcript_id IN ({}) 
    AND protein_variant LIKE '%F508%'
    LIMIT 10
    """.format(','.join(['?' for _ in test_transcripts]))
    
    result = pd.read_sql_query(query, conn, params=test_transcripts)
    print(f"\nFound {len(result)} F508 variants in CFTR transcripts:")
    print(result)
    
    # Get database statistics
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM alphamissense")
    total_rows = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT transcript_id) FROM alphamissense")
    unique_transcripts = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT protein_variant) FROM alphamissense")
    unique_variants = cursor.fetchone()[0]
    
    print(f"\nDatabase Statistics:")
    print(f"Total rows: {total_rows:,}")
    print(f"Unique transcripts: {unique_transcripts:,}")
    print(f"Unique protein variants: {unique_variants:,}")
    
    conn.close()
    print(f"\nDatabase created successfully: {db_file}")
    
    return db_file

def search_alphamissense(transcript_ids, protein_variant_pattern=None, limit=100):
    """
    Search AlphaMissense database for specific transcripts and optionally protein variants
    
    Args:
        transcript_ids (list): List of transcript IDs to search for
        protein_variant_pattern (str): Optional pattern to match in protein_variant column
        limit (int): Maximum number of results to return
    
    Returns:
        pandas.DataFrame: Search results
    """
    db_file = "cache/alphamissense/alphamissense_hg38.db"
    
    if not os.path.exists(db_file):
        print(f"Database file not found: {db_file}")
        print("Please run create_alphamissense_database() first")
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

if __name__ == "__main__":
    # Create the database
    db_file = create_alphamissense_database()
    
    # Example usage - search for F508 in CFTR transcripts
    print("\n" + "="*50)
    print("EXAMPLE SEARCH: F508 in CFTR transcripts")
    print("="*50)
    
    cftr_transcripts = [
        "ENST00000003084", "ENST00000426809", "ENST00000441019", "ENST00000472848",
        "ENST00000647720", "ENST00000647978", "ENST00000648260", "ENST00000649406",
        "ENST00000649781", "ENST00000685018", "ENST00000687278", "ENST00000699585",
        "ENST00000699596", "ENST00000699597", "ENST00000699598", "ENST00000699599",
        "ENST00000699600", "ENST00000699601", "ENST00000699602", "ENST00000699604",
        "ENST00000699605"
    ]
    
    results = search_alphamissense(cftr_transcripts, protein_variant_pattern="%F508%")
    if results is not None:
        print(f"Found {len(results)} F508 variants:")
        print(results) 