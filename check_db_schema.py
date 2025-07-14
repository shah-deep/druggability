#!/usr/bin/env python3
import sqlite3

conn = sqlite3.connect('cache/alphamissense/alphamissense_hg38.db')
cursor = conn.cursor()

# Get table info
cursor.execute('PRAGMA table_info(alphamissense)')
columns = cursor.fetchall()

print("AlphaMissense database columns:")
for col in columns:
    print(f"  {col[1]} ({col[2]})")

# Get a sample row to see the data structure
cursor.execute('SELECT * FROM alphamissense LIMIT 1')
sample = cursor.fetchone()
if sample:
    print(f"\nSample row: {sample}")

conn.close() 