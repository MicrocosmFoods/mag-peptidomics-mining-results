import pdfplumber
import csv
import argparse

def convert_pdf_table_to_tsv(pdf_path, output_path):
    try:
        # Open the PDF
        with pdfplumber.open(pdf_path) as pdf:
            # Open output TSV file
            with open(output_path, 'w', newline='') as tsv_file:
                writer = csv.writer(tsv_file, delimiter='\t')
                
                # Write header only once
                writer.writerow(['SAMPLE', 'PEPTIDE SEQUENCE', 'PROTEIN ACCESSION', 'DESCRIPTION'])
                
                # Initialize data_started flag outside the page loop
                data_started = False
                
                # Process all pages
                for page_num, page in enumerate(pdf.pages, 1):
                    print(f"Processing page {page_num}...")
                    # Extract text from current page
                    text = page.extract_text()
                    
                    if not text:
                        print(f"No text found on page {page_num}")
                        continue
                    
                    # Split into lines
                    lines = text.split('\n')
                    
                    # Process each line
                    for line in lines:
                        # Skip empty lines
                        if not line.strip():
                            continue
                        
                        # Check if we've reached the data section
                        if 'SAMPLE' in line and 'PEPTIDE SEQUENCE' in line:
                            data_started = True
                            continue
                        
                        # Process data lines
                        if data_started:
                            # Split the line into columns
                            parts = line.strip().split()
                            if len(parts) >= 4:
                                sample = parts[0]
                                description = parts[-1]
                                accession = parts[-2]
                                # Everything between sample and accession is the peptide sequence
                                peptide = ' '.join(parts[1:-2])
                                writer.writerow([sample, peptide, accession, description])
                    
                    print(f"Completed page {page_num}")
            
            print(f"Successfully converted table to {output_path}")
            print(f"Processed {len(pdf.pages)} pages")
            
    except Exception as e:
        print(f"Error processing PDF: {str(e)}")

def main():
    parser = argparse.ArgumentParser(
        description='Convert PDF table to TSV format, extracting SAMPLE, PEPTIDE SEQUENCE, PROTEIN ACCESSION, and DESCRIPTION columns from all pages.'
    )
    parser.add_argument(
        'input_pdf',
        help='Input PDF file containing the table'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output TSV file (default: input_basename_table.tsv)',
        default=None
    )
    
    args = parser.parse_args()
    
    # If no output file specified, create one based on input filename
    if args.output is None:
        output_path = args.input_pdf.split('.')[0] + '_table.tsv'
    else:
        output_path = args.output
        
    convert_pdf_table_to_tsv(args.input_pdf, output_path)

if __name__ == "__main__":
    main()