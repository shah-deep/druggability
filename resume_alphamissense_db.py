import pandas as pd
import sqlite3
import os
from pathlib import Path

def check_database_status():
    """
    Check the current status of the AlphaMissense database
    """
    db_file = "cache/alphamissense/alphamissense_hg38.db"
    
    if not os.path.exists(db_file):
        print("Database file does not exist yet.")
        return None, 0
    
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # Check if table exists
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='alphamissense'")
    if not cursor.fetchone():
        print("Database exists but table 'alphamissense' not found.")
        conn.close()
        return None, 0
    
    # Get current row count
    cursor.execute("SELECT COUNT(*) FROM alphamissense")
    current_rows = cursor.fetchone()[0]
    
    # Check if indexes exist
    cursor.execute("SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'")
    indexes = [row[0] for row in cursor.fetchall()]
    
    conn.close()
    
    print(f"Current database status:")
    print(f"  - Total rows: {current_rows:,}")
    print(f"  - Indexes: {indexes}")
    
    return db_file, current_rows

def count_tsv_rows():
    """
    Count total rows in the TSV file (excluding header comments)
    """
    tsv_file = "cache/alphamissense/AlphaMissense_hg38.tsv"
    
    if not os.path.exists(tsv_file):
        print(f"TSV file not found: {tsv_file}")
        return 0
    
    print("Counting rows in TSV file...")
    count = 0
    with open(tsv_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    
    print(f"Total data rows in TSV: {count:,}")
    return count

def resume_alphamissense_database():
    """
    Resume or create AlphaMissense database from TSV file
    """
    
    # File paths
    tsv_file = "cache/alphamissense/AlphaMissense_hg38.tsv"
    db_file = "cache/alphamissense/alphamissense_hg38.db"
    
    # Check current status
    existing_db, current_rows = check_database_status()
    total_tsv_rows = count_tsv_rows()
    
    if total_tsv_rows == 0:
        print("TSV file not found or empty. Cannot proceed.")
        return None
    
    # Calculate progress
    if current_rows > 0:
        progress = (current_rows / total_tsv_rows) * 100
        print(f"Progress: {progress:.1f}% ({current_rows:,}/{total_tsv_rows:,} rows)")
        
        if current_rows >= total_tsv_rows:
            print("Database appears to be complete!")
            return db_file
    else:
        print("Starting fresh database creation...")
    
    # Create database connection
    conn = sqlite3.connect(db_file)
    
    # Column names
    column_names = ['CHROM', 'POS', 'REF', 'ALT', 'genome', 'uniprot_id', 
                   'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class']
    
    # Process file in chunks
    chunk_size = 100000
    processed_rows = current_rows
    chunk_num = current_rows // chunk_size
    
    print(f"Processing file in chunks of {chunk_size} rows...")
    print(f"Starting from chunk {chunk_num + 1} (row {current_rows + 1})")
    
    # Skip already processed rows
    skip_rows = current_rows
    
    for chunk in pd.read_csv(tsv_file, 
                            sep='\t', 
                            comment='#',
                            names=column_names,
                            chunksize=chunk_size,
                            skiprows=skip_rows):
        
        # Clean up transcript_id column (remove version numbers)
        chunk['transcript_id'] = chunk['transcript_id'].str.split('.').str[0]
        
        # Write chunk to database
        chunk.to_sql('alphamissense', conn, if_exists='append', index=False)
        
        processed_rows += len(chunk)
        chunk_num += 1
        
        progress = (processed_rows / total_tsv_rows) * 100
        print(f"Processed chunk {chunk_num}: {len(chunk)} rows (Total: {processed_rows:,}, Progress: {progress:.1f}%)")
        
        # Show sample data from first chunk of this session
        if chunk_num == (current_rows // chunk_size) + 1:
            print(f"\nSample data from current chunk:")
            print(chunk.head())
    
    # Create indexes if they don't exist
    print("\nCreating/checking indexes...")
    
    indexes_to_create = [
        ("idx_transcript_id", "transcript_id"),
        ("idx_protein_variant", "protein_variant"),
        ("idx_uniprot_id", "uniprot_id")
    ]
    
    for index_name, column in indexes_to_create:
        try:
            conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON alphamissense({column})")
            print(f"Created/verified index: {index_name}")
        except Exception as e:
            print(f"Error creating index {index_name}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Final statistics
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM alphamissense")
    final_rows = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT transcript_id) FROM alphamissense")
    unique_transcripts = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT protein_variant) FROM alphamissense")
    unique_variants = cursor.fetchone()[0]
    
    print(f"\nFinal Database Statistics:")
    print(f"Total rows: {final_rows:,}")
    print(f"Unique transcripts: {unique_transcripts:,}")
    print(f"Unique protein variants: {unique_variants:,}")
    
    conn.close()
    print(f"\nDatabase completed successfully: {db_file}")
    
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
    results = search_alphamissense(cftr_transcripts, protein_variant_pattern="%F508%")
    if results is not None and len(results) > 0:
        print(f"Found {len(results)} F508 variants:")
        print(results)
    else:
        print("No F508 variants found in CFTR transcripts")

if __name__ == "__main__":
    # Check current status
    print("Checking current database status...")
    check_database_status()
    count_tsv_rows()
    
    # Resume or create database
    db_file = resume_alphamissense_database()
    
    if db_file:
        # Test the search function
        test_search_function() 