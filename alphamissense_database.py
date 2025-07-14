#!/usr/bin/env python3
"""
AlphaMissense Database Tool

This script provides efficient querying of AlphaMissense data by creating a database
from the large TSV file and enabling fast searches for specific transcripts and variants.
"""

import sqlite3
import pandas as pd
import json
import os
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AlphaMissenseDatabase:
    """Database class for efficient AlphaMissense data querying."""
    
    def __init__(self, tsv_file_path: str, db_path: str = "alphamissense.db"):
        """
        Initialize the AlphaMissense database.
        
        Args:
            tsv_file_path: Path to the AlphaMissense TSV file
            db_path: Path for the SQLite database file
        """
        self.tsv_file_path = tsv_file_path
        self.db_path = db_path
        self.conn = None
        
    def create_database(self, chunk_size: int = 10000) -> None:
        """
        Create SQLite database from AlphaMissense TSV file.
        
        Args:
            chunk_size: Number of rows to process at once (for memory efficiency)
        """
        logger.info(f"Creating database from {self.tsv_file_path}")
        
        # Connect to SQLite database
        self.conn = sqlite3.connect(self.db_path)
        
        # Read TSV file in chunks and insert into database
        chunk_count = 0
        for chunk in pd.read_csv(self.tsv_file_path, sep='\t', chunksize=chunk_size):
            chunk_count += 1
            logger.info(f"Processing chunk {chunk_count}")
            
            # Create table if first chunk
            if chunk_count == 1:
                chunk.to_sql('alphamissense', self.conn, if_exists='replace', index=False)
            else:
                chunk.to_sql('alphamissense', self.conn, if_exists='append', index=False)
        
        # Create indexes for faster querying
        logger.info("Creating indexes for faster querying...")
        self.conn.execute("CREATE INDEX IF NOT EXISTS idx_transcript ON alphamissense(transcript)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS idx_position ON alphamissense(position)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS idx_gene ON alphamissense(gene)")
        
        logger.info(f"Database created successfully: {self.db_path}")
        
    def query_transcripts(self, transcript_ids: List[str], additional_filters: Optional[Dict] = None) -> pd.DataFrame:
        """
        Query AlphaMissense data for specific transcript IDs.
        
        Args:
            transcript_ids: List of transcript IDs to search for
            additional_filters: Additional filters to apply (e.g., {'position': 'F508'})
            
        Returns:
            DataFrame with matching results
        """
        if not self.conn:
            raise ValueError("Database not initialized. Call create_database() first.")
        
        # Build query
        transcript_list = "', '".join(transcript_ids)
        query = f"SELECT * FROM alphamissense WHERE transcript IN ('{transcript_list}')"
        
        # Add additional filters
        if additional_filters:
            for key, value in additional_filters.items():
                if isinstance(value, str):
                    query += f" AND {key} LIKE '%{value}%'"
                else:
                    query += f" AND {key} = '{value}'"
        
        logger.info(f"Executing query: {query}")
        return pd.read_sql_query(query, self.conn)
    
    def search_by_pattern(self, pattern: str, column: str = "transcript") -> pd.DataFrame:
        """
        Search AlphaMissense data using a pattern.
        
        Args:
            pattern: Pattern to search for
            column: Column to search in
            
        Returns:
            DataFrame with matching results
        """
        if not self.conn:
            raise ValueError("Database not initialized. Call create_database() first.")
        
        query = f"SELECT * FROM alphamissense WHERE {column} LIKE '%{pattern}%'"
        logger.info(f"Executing pattern search: {query}")
        return pd.read_sql_query(query, self.conn)
    
    def get_database_info(self) -> Dict:
        """Get information about the database."""
        if not self.conn:
            return {"error": "Database not initialized"}
        
        cursor = self.conn.cursor()
        
        # Get table info
        cursor.execute("PRAGMA table_info(alphamissense)")
        columns = cursor.fetchall()
        
        # Get row count
        cursor.execute("SELECT COUNT(*) FROM alphamissense")
        row_count = cursor.fetchone()[0]
        
        return {
            "database_path": self.db_path,
            "table_name": "alphamissense",
            "row_count": row_count,
            "columns": [col[1] for col in columns]
        }
    
    def close(self):
        """Close the database connection."""
        if self.conn:
            self.conn.close()
            logger.info("Database connection closed")

def load_processed_input(input_file: str) -> Dict:
    """Load processed input JSON file."""
    with open(input_file, 'r') as f:
        return json.load(f)

def extract_transcript_ids(data: Dict) -> List[str]:
    """Extract all transcript IDs from processed input data."""
    transcript_ids = set()
    
    for variant in data.get('missense_variants', []):
        # Add main transcript ID
        if variant.get('transcript_id'):
            transcript_ids.add(variant['transcript_id'])
        
        # Add VEP transcript IDs
        for transcript_id in variant.get('vep_transcript_ids', []):
            transcript_ids.add(transcript_id)
    
    return list(transcript_ids)

def main():
    """Main function to demonstrate usage."""
    
    # Configuration
    tsv_file_path = "cache/alphamissense/AlphaMissense_hg38.tsv"
    db_path = "alphamissense.db"
    input_file = "processed_input.json"
    
    # Check if TSV file exists
    if not os.path.exists(tsv_file_path):
        logger.error(f"AlphaMissense TSV file not found at: {tsv_file_path}")
        logger.info("Please ensure the file exists or update the path in the script.")
        return
    
    # Initialize database
    db = AlphaMissenseDatabase(tsv_file_path, db_path)
    
    # Create database if it doesn't exist
    if not os.path.exists(db_path):
        logger.info("Creating new database...")
        db.create_database()
    else:
        logger.info("Using existing database")
        db.conn = sqlite3.connect(db_path)
    
    # Load processed input data
    try:
        input_data = load_processed_input(input_file)
        transcript_ids = extract_transcript_ids(input_data)
        logger.info(f"Found {len(transcript_ids)} transcript IDs in processed input")
        
        # Example 1: Query specific transcripts
        logger.info("Querying specific transcripts...")
        results = db.query_transcripts(transcript_ids[:5])  # First 5 for demo
        print(f"Found {len(results)} records for specified transcripts")
        if not results.empty:
            print(results.head())
        
        # Example 2: Search for F508 pattern (CFTR variant)
        logger.info("Searching for F508 pattern...")
        f508_results = db.search_by_pattern("F508")
        print(f"Found {len(f508_results)} records containing F508")
        if not f508_results.empty:
            print(f508_results.head())
        
        # Example 3: Combined search (transcripts + F508)
        logger.info("Combined search: transcripts + F508 pattern...")
        cftr_transcripts = [tid for tid in transcript_ids if 'ENST00000003084' in tid or 'CFTR' in tid]
        combined_results = db.query_transcripts(cftr_transcripts, {'position': 'F508'})
        print(f"Found {len(combined_results)} records for CFTR transcripts with F508")
        if not combined_results.empty:
            print(combined_results.head())
        
        # Get database info
        info = db.get_database_info()
        print(f"\nDatabase Info: {info}")
        
    except Exception as e:
        logger.error(f"Error processing data: {e}")
    
    finally:
        db.close()

if __name__ == "__main__":
    main() 