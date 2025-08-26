import argparse
import os

def parse_pride_tsv(input_file):
    # Get output filename by removing extension and adding '_parsed.tsv'
    filename = os.path.basename(input_file)
    base_name = filename.split('.')[0]
    output_file = f"{base_name}_parsed.tsv"
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Find the PSH header line and write it
        for line in infile:
            if line.startswith('PSH'):
                outfile.write(line)
                break
        
        # Write all PSM data lines
        for line in infile:
            if line.startswith('PSM'):
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(
        description='Parse PRIDE TSV files to extract PSH header and PSM data lines.'
    )
    parser.add_argument(
        'input_file',
        help='Input PRIDE TSV file to parse'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file name (default: input_basename_parsed.tsv)',
        default=None
    )

    args = parser.parse_args()
    parse_pride_tsv(args.input_file)

if __name__ == "__main__":
    main()