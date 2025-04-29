#!/usr/bin/env python3
"""
Combine multiple batches of mining results into single files.

This script takes multiple directories of batch runs and creates combined results.
For TSV files, it adds a 'batch' column to track the source batch.
For FASTA files, it validates the format and combines sequences.
Only all_molecule_counts.tsv will be checked for duplicates by genome_name.
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
import re
from collections import defaultdict

def extract_batch_name(directory):
    """Extract the batch name from the directory path."""
    # Extract everything before 'main-results'
    match = re.search(r'(.+?)-main-results', directory)
    if match:
        return match.group(1).split('/')[-1]  # Get the last part of the path
    return os.path.basename(directory)  # Fallback to directory name

def is_fasta_file(file_path):
    """Check if a file is in FASTA format."""
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        return first_line.startswith('>')

def combine_tsv_files(file_paths, output_file, batch_names, check_duplicates=False):
    """
    Combine multiple TSV files into one.
    
    Args:
        file_paths: List of paths to TSV files
        output_file: Path to save the combined file
        batch_names: List of batch names corresponding to file_paths
        check_duplicates: Whether to check for duplicates in the first column
    """
    if not file_paths:
        print(f"No files to combine for {output_file}")
        return
    
    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Read all dataframes
    dfs = []
    for i, file_path in enumerate(file_paths):
        try:
            df = pd.read_csv(file_path, sep='\t')
            # Add batch column if it doesn't exist
            if 'batch' not in df.columns:
                df['batch'] = batch_names[i]
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    
    if not dfs:
        print(f"No valid data found for {output_file}")
        return
    
    # Check column consistency
    columns = dfs[0].columns
    for i, df in enumerate(dfs[1:], 1):
        if set(df.columns) != set(columns):
            print(f"Warning: Column mismatch in {file_paths[i]}")
            print(f"  Expected: {columns}")
            print(f"  Found: {df.columns}")
    
    # Combine dataframes
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Check for duplicates if needed
    if check_duplicates and len(combined_df) > 0:
        id_col = combined_df.columns[0]  # First column is the ID column
        duplicates = combined_df[id_col].duplicated()
        if duplicates.any():
            print(f"Found {duplicates.sum()} duplicates in {output_file}")
            print(f"Keeping first occurrence of each {id_col}")
            combined_df = combined_df.drop_duplicates(subset=[id_col], keep='first')
    
    # Save combined dataframe
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Combined {len(file_paths)} files into {output_file}")

def combine_fasta_files(file_paths, output_file):
    """
    Combine multiple FASTA files into one, checking for duplicate headers.
    
    Args:
        file_paths: List of paths to FASTA files
        output_file: Path to save the combined file
    """
    if not file_paths:
        print(f"No files to combine for {output_file}")
        return
    
    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Track headers to check for duplicates
    seen_headers = set()
    duplicate_count = 0
    
    with open(output_file, 'w') as out_f:
        for file_path in file_paths:
            try:
                with open(file_path, 'r') as in_f:
                    current_header = None
                    for line in in_f:
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('>'):
                            current_header = line
                            if current_header in seen_headers:
                                duplicate_count += 1
                                # Skip this sequence
                                current_header = None
                            else:
                                seen_headers.add(current_header)
                                out_f.write(line + '\n')
                        elif current_header is not None:
                            out_f.write(line + '\n')
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
    
    if duplicate_count > 0:
        print(f"Found and skipped {duplicate_count} duplicate headers in {output_file}")
    
    print(f"Combined {len(file_paths)} files into {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Combine batch mining results')
    parser.add_argument('input_dirs', nargs='+', help='Input directories containing batch results')
    parser.add_argument('--output-dir', required=True, help='Output directory for combined results')
    args = parser.parse_args()
    
    # Validate input directories
    valid_dirs = []
    batch_names = []
    for dir_path in args.input_dirs:
        if not os.path.isdir(dir_path):
            print(f"Warning: {dir_path} is not a directory, skipping")
            continue
        valid_dirs.append(dir_path)
        batch_names.append(extract_batch_name(dir_path))
    
    if not valid_dirs:
        print("Error: No valid input directories provided")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Collect files to combine
    tsv_files = defaultdict(list)
    fasta_files = defaultdict(list)
    
    # Process main directory files
    for dir_path in valid_dirs:
        for file_name in os.listdir(dir_path):
            file_path = os.path.join(dir_path, file_name)
            if os.path.isfile(file_path) and not os.path.isdir(file_path):
                if file_name.endswith('.tsv'):
                    tsv_files[file_name].append(file_path)
                elif is_fasta_file(file_path):
                    fasta_files[file_name].append(file_path)
    
    # Process bgc_info subdirectory files
    for dir_path in valid_dirs:
        bgc_info_dir = os.path.join(dir_path, 'bgc_info')
        if os.path.isdir(bgc_info_dir):
            for file_name in os.listdir(bgc_info_dir):
                file_path = os.path.join(bgc_info_dir, file_name)
                if os.path.isfile(file_path):
                    if file_name.endswith('.tsv'):
                        tsv_files[os.path.join('bgc_info', file_name)].append(file_path)
                    elif is_fasta_file(file_path):
                        fasta_files[os.path.join('bgc_info', file_name)].append(file_path)
    
    # Combine TSV files
    for file_name, file_paths in tsv_files.items():
        output_file = os.path.join(args.output_dir, file_name)
        batch_indices = [valid_dirs.index(os.path.dirname(file_path.replace('/bgc_info', ''))) 
                         for file_path in file_paths]
        file_batch_names = [batch_names[i] for i in batch_indices]
        
        # Only check for duplicates in all_molecule_counts.tsv
        check_duplicates = (file_name == 'all_molecule_counts.tsv')
        combine_tsv_files(file_paths, output_file, file_batch_names, check_duplicates)
    
    # Combine FASTA files
    for file_name, file_paths in fasta_files.items():
        output_file = os.path.join(args.output_dir, file_name)
        combine_fasta_files(file_paths, output_file)
    
    print(f"Combined results saved to {args.output_dir}")

if __name__ == "__main__":
    main()
