#!/usr/bin/env python3
"""
Database Indexing and Optimization Script

This script examines the current ClinVar cache database structure,
removes unnecessary indexes, and keeps only the primary key index.
"""

import sqlite3
import logging
from pathlib import Path
from typing import List, Dict, Any

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def examine_database_structure(db_path: str) -> Dict[str, Any]:
    """Examine the current database structure and indexes"""
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Get table information
        cursor.execute("PRAGMA table_info(variant_annotations)")
        columns = cursor.fetchall()
        
        # Get current indexes
        cursor.execute("PRAGMA index_list(variant_annotations)")
        indexes = cursor.fetchall()
        
        # Get table statistics
        cursor.execute("SELECT COUNT(*) FROM variant_annotations")
        row_count = cursor.fetchone()[0]
        
        conn.close()
        
        return {
            'columns': columns,
            'indexes': indexes,
            'row_count': row_count
        }
        
    except Exception as e:
        logger.error(f"Error examining database structure: {e}")
        return {}


def remove_unnecessary_indexes(db_path: str) -> bool:
    """Remove all indexes except the primary key index"""
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Enable foreign keys and WAL mode for better performance
        cursor.execute("PRAGMA foreign_keys = ON")
        cursor.execute("PRAGMA journal_mode = WAL")
        cursor.execute("PRAGMA synchronous = NORMAL")
        cursor.execute("PRAGMA cache_size = 10000")
        cursor.execute("PRAGMA temp_store = MEMORY")
        
        # Get all current indexes
        cursor.execute("PRAGMA index_list(variant_annotations)")
        indexes = cursor.fetchall()
        
        # Remove all indexes except the primary key auto-index
        for index_info in indexes:
            index_name = index_info[1]
            # Keep only the primary key auto-index
            if not index_name.startswith('sqlite_autoindex'):
                logger.info(f"Dropping index: {index_name}")
                cursor.execute(f"DROP INDEX IF EXISTS {index_name}")
        
        # Commit changes
        conn.commit()
        conn.close()
        
        logger.info("Successfully removed unnecessary indexes")
        return True
        
    except Exception as e:
        logger.error(f"Error removing indexes: {e}")
        return False


def test_database_performance(db_path: str) -> Dict[str, Any]:
    """Test database performance with primary key lookups"""
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        results = {}
        
        # Test 1: Primary key lookup (this is the main lookup pattern)
        cursor.execute("SELECT COUNT(*) FROM variant_annotations WHERE variant_key LIKE '%TP53%'")
        results['primary_key_lookup'] = cursor.fetchone()[0]
        
        # Test 2: Query execution time for primary key lookup
        import time
        start_time = time.time()
        cursor.execute("SELECT * FROM variant_annotations WHERE variant_key = 'test_key'")
        cursor.fetchall()
        results['primary_key_query_time_ms'] = (time.time() - start_time) * 1000
        
        # Test 3: Full table scan (for comparison)
        start_time = time.time()
        cursor.execute("SELECT COUNT(*) FROM variant_annotations")
        cursor.fetchall()
        results['full_table_scan_time_ms'] = (time.time() - start_time) * 1000
        
        conn.close()
        
        return results
        
    except Exception as e:
        logger.error(f"Error testing database performance: {e}")
        return {}


def analyze_query_plans(db_path: str) -> List[str]:
    """Analyze query execution plans to verify primary key index usage"""
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        plans = []
        
        # Test queries with EXPLAIN QUERY PLAN
        test_queries = [
            "SELECT * FROM variant_annotations WHERE variant_key = 'test_key'",
            "SELECT * FROM variant_annotations WHERE variant_key LIKE '%TP53%'",
            "SELECT * FROM variant_annotations WHERE variant_key = 'TP53_p.Arg175His_17_7675088_G_A'"
        ]
        
        for i, query in enumerate(test_queries, 1):
            cursor.execute(f"EXPLAIN QUERY PLAN {query}")
            plan = cursor.fetchall()
            plans.append(f"Query {i}: {query}")
            plans.append("Plan:")
            for line in plan:
                plans.append(f"  {line[3]}")
            plans.append("")
        
        conn.close()
        return plans
        
    except Exception as e:
        logger.error(f"Error analyzing query plans: {e}")
        return []


def main():
    """Main function to examine, optimize, and test the database"""
    db_path = "cache/clinvar/clinvar_cache.db"
    
    print("=" * 80)
    print("DATABASE INDEX OPTIMIZATION - PRIMARY KEY ONLY")
    print("=" * 80)
    
    # Step 1: Examine current structure
    print("\n1. EXAMINING CURRENT DATABASE STRUCTURE")
    print("-" * 50)
    
    structure = examine_database_structure(db_path)
    if structure:
        print(f"Database has {structure['row_count']} rows")
        print(f"Number of columns: {len(structure['columns'])}")
        print(f"Number of existing indexes: {len(structure['indexes'])}")
        
        print("\nColumns:")
        for col in structure['columns']:
            print(f"  - {col[1]} ({col[2]})")
        
        print("\nExisting indexes:")
        for idx in structure['indexes']:
            print(f"  - {idx[1]}")
    else:
        print("Could not examine database structure")
    
    # Step 2: Remove unnecessary indexes
    print("\n2. REMOVING UNNECESSARY INDEXES")
    print("-" * 50)
    
    if remove_unnecessary_indexes(db_path):
        print("✓ Successfully removed unnecessary indexes")
    else:
        print("✗ Failed to remove indexes")
    
    # Step 3: Re-examine structure after removing indexes
    print("\n3. VERIFYING INDEX REMOVAL")
    print("-" * 50)
    
    structure_after = examine_database_structure(db_path)
    if structure_after:
        print(f"Number of indexes after optimization: {len(structure_after['indexes'])}")
        print("\nRemaining indexes:")
        for idx in structure_after['indexes']:
            print(f"  - {idx[1]}")
    
    # Step 4: Test performance
    print("\n4. TESTING DATABASE PERFORMANCE")
    print("-" * 50)
    
    performance_results = test_database_performance(db_path)
    if performance_results:
        print("Performance test results:")
        for test, result in performance_results.items():
            print(f"  {test}: {result}")
    
    # Step 5: Analyze query plans
    print("\n5. ANALYZING QUERY EXECUTION PLANS")
    print("-" * 50)
    
    query_plans = analyze_query_plans(db_path)
    if query_plans:
        for plan in query_plans:
            print(plan)
    
    print("\n" + "=" * 80)
    print("DATABASE OPTIMIZATION COMPLETE - PRIMARY KEY ONLY")
    print("=" * 80)


if __name__ == "__main__":
    main() 