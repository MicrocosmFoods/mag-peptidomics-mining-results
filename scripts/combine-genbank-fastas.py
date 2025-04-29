#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_genome_scaffold_map(tsv_file):
    """
    Load genome to scaffold mapping from TSV file.
    
    Args:
        tsv_file (Path): Path to TSV file with genome and scaffold names
        
    Returns:
        dict: Mapping of scaffold names to genome names
    """
    df = pd.read_csv(tsv_file, sep='\t', header=None)
    return dict(zip(df[1], df[0]))

def clean_and_combine_fastas(input_files, output_file, genome_tsv=None):
    """
    Combine multiple FASTA files and clean headers to keep only text before first whitespace.
    If genome_tsv is provided, fix antismash headers using genome-scaffold mapping.
    
    Args:
        input_files (list): List of input FASTA file paths
        output_file (str): Path to output combined FASTA
        genome_tsv (Path, optional): Path to genome-scaffold mapping TSV
    """
    # Load genome-scaffold mapping if provided
    scaffold_to_genome = None
    if genome_tsv:
        scaffold_to_genome = load_genome_scaffold_map(genome_tsv)
    
    # Store all cleaned records
    combined_records = []
    
    # Process each input file
    for fasta_file in input_files:
        # Open and parse the file
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Clean header - keep only text before first whitespace
                clean_id = record.id.split()[0]
                
                # If this is the antismash file and we have genome mapping
                if scaffold_to_genome and fasta_file.name == "antismash_peptides.fasta":
                    scaffold_name = clean_id.split('_')[0]  # Get scaffold name
                    if scaffold_name in scaffold_to_genome:
                        genome_name = scaffold_to_genome[scaffold_name]
                        clean_id = f"{genome_name}_id_{clean_id}"
                
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
    
    parser.add_argument(
        '-g', '--genome_tsv',
        help='TSV file mapping genome names to scaffold names (optional)'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert input files to Path objects and check they exist
    input_files = [Path(f) for f in args.input]
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file not found: {f}")
    
    # Check genome TSV if provided
    genome_tsv = None
    if args.genome_tsv:
        genome_tsv = Path(args.genome_tsv)
        if not genome_tsv.exists():
            raise FileNotFoundError(f"Genome TSV not found: {genome_tsv}")
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Process files
    clean_and_combine_fastas(input_files, output_path, genome_tsv)

if __name__ == "__main__":
    main()