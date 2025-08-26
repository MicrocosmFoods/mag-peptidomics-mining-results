import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd

def pride_to_fasta(input_tsv, output_fasta, sample_name=None, seq_column='sequence'):
    # Read TSV file
    try:
        df = pd.read_csv(input_tsv, sep='\t')
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return
    
    # Check if sequence column exists
    if seq_column not in df.columns:
        raise ValueError(f"Column '{seq_column}' not found in TSV. Available columns: {', '.join(df.columns)}")
    
    # Get base filename without extension (remove all extensions after first '.')
    base_filename = input_tsv.split('.')[0]
    base_filename = os.path.basename(base_filename)
    
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
        description='Convert PRIDE TSV file with peptide sequences to FASTA format'
    )
    parser.add_argument(
        'input_tsv',
        help='Input PRIDE TSV file containing peptide sequences'
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
        help='Name of column containing peptide sequences (default: "sequence")',
        default='sequence'
    )
    
    args = parser.parse_args()
    
    # If no output file specified, create one based on input filename
    if args.output is None:
        output_path = args.input_tsv.split('.')[0] + '.fasta'
    else:
        output_path = args.output
        
    pride_to_fasta(args.input_tsv, output_path, args.sample, args.column)

if __name__ == "__main__":
    main()