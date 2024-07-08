import os
import argparse
import logging
from Bio import Entrez

def download_fasta_gff(accession, email, output_dir):
    try:
        Entrez.email = email

        # Fetch the FASTA file
        logging.info(f"Downloading FASTA file for {accession}...")
        fasta_handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
        fasta_data = fasta_handle.read()
        fasta_handle.close()

        # Save the FASTA file
        fasta_file = os.path.join(output_dir, f"{accession}.fasta")
        with open(fasta_file, 'w') as f:
            f.write(fasta_data)
        logging.info(f"FASTA file saved as {fasta_file}")

        # Fetch the GFF file
        logging.info(f"Downloading GFF file for {accession}...")
        gff_handle = Entrez.efetch(db="nuccore", id=accession, rettype="gff3", retmode="text")
        gff_data = gff_handle.read()
        gff_handle.close()

        # Save the GFF file
        gff_file = os.path.join(output_dir, f"{accession}.gff")
        with open(gff_file, 'w') as f:
            f.write(gff_data)
        logging.info(f"GFF file saved as {gff_file}")

        return fasta_file, gff_file

    except Exception as e:
        logging.error(f"Error downloading files for {accession}: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download FASTA and GFF files from NCBI.")
    parser.add_argument("-a", "--accession", required=True, help="NCBI accession number")
    parser.add_argument("-e", "--email", default="your_default_email@example.com", help="Email address (required by NCBI)")
    parser.add_argument("-o", "--outputdir", default=".", help="Output directory (default is current directory)")

    args = parser.parse_args()
    
    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)
    
    logging.basicConfig(filename=f'{args.accession}_DataRetrieval.log', level=logging.INFO)
    
    download_fasta_gff(args.accession, args.email, args.outputdir)
