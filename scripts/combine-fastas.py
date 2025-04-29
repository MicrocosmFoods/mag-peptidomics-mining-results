#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def clean_and_combine_fastas(input_files, output_file):
    """
    Combine multiple FASTA files and clean headers to keep only text before first whitespace.
    
    Args:
        input_files (list): List of input FASTA file paths
        output_file (str): Path to output combined FASTA
    """
    # Store all cleaned records
    combined_records = []
    
    # Process each input file
    for fasta_file in input_files:
        # Open and parse the file
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Clean header - keep only text before first whitespace
                clean_id = record.id.split()[0]
                
                # Create new record with cleaned header
                new_record = SeqRecord(
                    seq=record.seq,
                    id=clean_id,
                    description=""  # Empty description to avoid redundant header info
                )
                
                combined_records.append(new_record)
    
    # Write combined records to output file
    with open(output_file, "w") as out_handle:
        SeqIO.write(combined_records, out_handle, "fasta")
    
    # Print summary
    print(f"Combined {len(input_files)} files into {output_file}")
    print(f"Total sequences: {len(combined_records)}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Combine multiple FASTA files and clean headers'
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        nargs='+',
        help='Input FASTA files (space-separated)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output combined FASTA file'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert input files to Path objects and check they exist
    input_files = [Path(f) for f in args.input]
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file not found: {f}")
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Process files
    clean_and_combine_fastas(input_files, output_path)

if __name__ == "__main__":
    main()