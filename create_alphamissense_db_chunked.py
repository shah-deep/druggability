import pandas as pd
import sqlite3
import os
from pathlib import Path

def create_alphamissense_database_chunked():
    """
    Convert AlphaMissense_hg38.tsv to SQLite database using chunked processing
    """
    
    # File paths
    tsv_file = "cache/alphamissense/AlphaMissense_hg38.tsv"
    db_file = "cache/alphamissense/alphamissense_hg38.db"
    
    print(f"Reading TSV file in chunks: {tsv_file}")
    
    # Create database connection
    conn = sqlite3.connect(db_file)
    
    # Column names based on the file structure
    column_names = ['CHROM', 'POS', 'REF', 'ALT', 'genome', 'uniprot_id', 
                   'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class']
    
    # Process file in chunks
    chunk_size = 100000  # Process 100k rows at a time
    total_rows = 0
    
    # First, let's check the file structure
    print("Checking file structure...")
    with open(tsv_file, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            if i == 0:  # First non-comment line
                print(f"First data line: {line.strip()}")
                break
    
    # Process the file in chunks
    print(f"Processing file in chunks of {chunk_size} rows...")
    
    for chunk_num, chunk in enumerate(pd.read_csv(tsv_file, 
                                                 sep='\t', 
                                                 comment='#',
                                                 names=column_names,
                                                 chunksize=chunk_size)):
        
        # Clean up transcript_id column (remove version numbers)
        chunk['transcript_id'] = chunk['transcript_id'].str.split('.').str[0]
        
        # Write chunk to database
        if chunk_num == 0:
            # First chunk - create table
            chunk.to_sql('alphamissense', conn, if_exists='replace', index=False)
        else:
            # Subsequent chunks - append to table
            chunk.to_sql('alphamissense', conn, if_exists='append', index=False)
        
        total_rows += len(chunk)
        print(f"Processed chunk {chunk_num + 1}: {len(chunk)} rows (Total: {total_rows:,})")
        
        # Show sample data from first chunk
        if chunk_num == 0:
            print(f"\nSample data from first chunk:")
            print(chunk.head())
    
    # Create indexes for fast searching
    print("\nCreating indexes...")
    
    # Index on transcript_id for fast transcript searches
    conn.execute("CREATE INDEX IF NOT EXISTS idx_transcript_id ON alphamissense(transcript_id)")
    print("Created transcript_id index")
    
    # Index on protein_variant for fast protein variant searches
    conn.execute("CREATE INDEX IF NOT EXISTS idx_protein_variant ON alphamissense(protein_variant)")
    print("Created protein_variant index")
    
    # Index on uniprot_id for protein searches
    conn.execute("CREATE INDEX IF NOT EXISTS idx_uniprot_id ON alphamissense(uniprot_id)")
    print("Created uniprot_id index")
    
    # Index on chromosome and position for genomic searches
    # conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos ON alphamissense(CHROM, POS)")
    # print("Created chrom_pos index")
    
    # # Index on pathogenicity class
    # conn.execute("CREATE INDEX IF NOT EXISTS idx_am_class ON alphamissense(am_class)")
    # print("Created am_class index")
    
    # Commit changes
    conn.commit()
    
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
        print("Please run create_alphamissense_database_chunked() first")
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
    results = search_alphamissense(cftr_transcripts, protein_variant_pattern="%F508%")
    if results is not None:
        print(f"Found {len(results)} F508 variants:")
        print(results)
    else:
        print("No results found or database not created yet")

if __name__ == "__main__":
    # Create the database
    db_file = create_alphamissense_database_chunked()
    
    # Test the search function
    test_search_function() 