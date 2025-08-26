import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd

def csv_to_fasta(input_csv, output_fasta, sample_name=None, seq_column='Peptide sequence'):
    # Read CSV file
    df = pd.read_csv(input_csv)
    
    # Check if sequence column exists
    if seq_column not in df.columns:
        raise ValueError(f"Column '{seq_column}' not found in CSV. Available columns: {', '.join(df.columns)}")
    
    # Get base filename without extension
    base_filename = os.path.splitext(os.path.basename(input_csv))[0]
    
    # Create sequence records
    records = []
    for idx, row in enumerate(df[seq_column], 1):
        # Create header ID
        if sample_name:
            seq_id = f"{sample_name}_{base_filename}_peptide_{idx}"
        else:
            seq_id = f"{base_filename}_peptide_{idx}"
            
        # Create sequence record
        sequence = Seq(row)
        record = SeqRecord(sequence,
                          id=seq_id,
                          description="")
        records.append(record)
    
    # Write to FASTA file
    SeqIO.write(records, output_fasta, "fasta")
    print(f"Converted {len(records)} sequences to {output_fasta}")

def main():
    parser = argparse.ArgumentParser(
        description='Convert CSV file with peptide sequences to FASTA format'
    )
    parser.add_argument(
        'input_csv',
        help='Input CSV file containing peptide sequences'
    )
    parser.add_argument(
        '-s', '--sample',
        help='Optional sample name to include in sequence headers',
        default=None
    )
    parser.add_argument(
        '-o', '--output',
        help='Output FASTA file (default: input_basename.fasta)',
        default=None
    )
    parser.add_argument(
        '-c', '--column',
        help='Name of column containing peptide sequences (default: "Peptide sequence")',
        default='Peptide sequence'
    )
    
    args = parser.parse_args()
    
    # If no output file specified, create one based on input filename
    if args.output is None:
        output_path = os.path.splitext(args.input_csv)[0] + '.fasta'
    else:
        output_path = args.output
        
    csv_to_fasta(args.input_csv, output_path, args.sample, args.column)

if __name__ == "__main__":
    main()