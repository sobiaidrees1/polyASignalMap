import os
import argparse
import logging
import subprocess
import shutil
import glob
import sys

# Define the default polyA signals
default_polyA_signals = ["AATAAA"]

def run_command(command):
    try:
        logging.info(f"Running command: {command}")
        subprocess.check_call(command, shell=True)
        logging.info("Command execution completed.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e.cmd}")
        logging.error(f"Return code: {e.returncode}")
        raise

def overwrite_directory(output_dir):
    if os.path.exists(output_dir):
        user_input = input(f"Directory '{output_dir}' already exists. Do you want to overwrite it? (y/n): ")
        if user_input.lower() != 'y':
            print("Please provide a different output directory.")
            sys.exit(1)
        logging.warning(f"Deleting existing directory '{output_dir}' and its contents to overwrite.")
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

def get_fastq_files(input_dir):
    # Get all FASTQ files in the input directory
    fastq_files = glob.glob(os.path.join(input_dir, "*.fq.gz")) + glob.glob(os.path.join(input_dir, "*.fastq.gz"))
    return fastq_files

def get_paired_files(input_dir, extension):
    # Get all R1 files
    r1_files = glob.glob(os.path.join(input_dir, "*_R1*.fq.gz")) + glob.glob(os.path.join(input_dir, "*_R1*.fastq.gz"))
    paired_files = []

    for r1 in r1_files:
        if extension == 'R2':
            r2 = r1.replace('_R1', '_R2')
        elif extension == 'R3':
            r2 = r1.replace('_R1', '_R3')
        if os.path.exists(r2):
            paired_files.append((r1, r2))
        else:
            logging.warning(f"Paired file not found for {r1}, expected {r2}")

    return paired_files

if __name__ == "__main__":
    default_email = "your_default_email@example.com"

    parser = argparse.ArgumentParser(description="Pipeline for downloading, processing polyA sequences, and aligning with STAR.")
    parser.add_argument("-a", "--accession", required=True, help="NCBI accession number")
    parser.add_argument("-o", "--outputdir", default=".", help="Output directory (default is current directory)")
    parser.add_argument("-m", "--max_read_length", type=int, default=100, help="Maximum read length for STAR indexing (default is 100)")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing FASTQ files")
    parser.add_argument("-u", "--upstream_length", type=int, default=20, help="Length of upstream sequence to extract (default is 20)")
    parser.add_argument("-t", "--type", choices=['paired', 'single'], help="Type of sequencing data (paired-end or single-end)")
    parser.add_argument("-p", "--pipeline_steps", nargs='+', choices=['download', 'convert', 'polyA', 'align'], default=['download', 'convert', 'polyA', 'align'], help="Pipeline steps to execute (default: all)")
    parser.add_argument("-x", "--extension", choices=['R2', 'R3'], default='R2', help="Extension for the second file in paired-end data (default: R2)")
    parser.add_argument("-s", "--polyA_signals", nargs='+', default=default_polyA_signals, help="PolyA signal sequence to detect (default: AATAAA)")
    parser.add_argument("-n", "--genomeSAindexNbases", type=int, default=6, help="Number of bases to be used for STAR genome indexing (default is 6)")

    args = parser.parse_args()

    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    logging.basicConfig(filename=f'{args.accession}_main.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    print(f"Searching for polyA signals: {', '.join(args.polyA_signals)}")

    try:
        if 'download' in args.pipeline_steps:
            logging.info("Step 1: Starting download of FASTA and GFF files.")
            run_command(f"python3 src/dataretrieval.py -a {args.accession} -o {args.outputdir}")

        if 'convert' in args.pipeline_steps:
            gff_file = os.path.join(args.outputdir, f"{args.accession}.gff")
            gtf_file = os.path.join(args.outputdir, f"{args.accession}.gtf")
            logging.info("Step 2: Starting conversion of GFF to GTF.")
            run_command(f"python3 src/coordinates.py -g {gff_file} -o {gtf_file}")

        if 'polyA' in args.pipeline_steps:
            fasta_file = os.path.join(args.outputdir, f"{args.accession}.fasta")
            gtf_file = os.path.join(args.outputdir, f"{args.accession}.gtf")
            logging.info("Step 3: Starting polyA sequence processing.")
            run_command(f"python3 src/processing.py -f {fasta_file} -g {gtf_file} -o {args.outputdir} -a {args.accession} -s {' '.join(args.polyA_signals)}")

            gtf_file = os.path.join(args.outputdir, f"{args.accession}_processed.gtf")

            indexed_genome_path = os.path.join(args.outputdir, "indexed_genome")
            if not os.path.exists(indexed_genome_path):
                os.makedirs(indexed_genome_path)

            genome_file = os.path.join(args.outputdir, f"{args.accession}.fasta")
            logging.info("Step 4: Starting STAR indexing.")
            command = (
                f"STAR "
                f"--runMode genomeGenerate "
                f"--genomeDir {indexed_genome_path} "
                f"--genomeFastaFiles {genome_file} "
                f"--sjdbGTFfile {gtf_file} "
                f"--sjdbOverhang {args.max_read_length - 1} "
                f"--genomeSAindexNbases {args.genomeSAindexNbases} "  # Use user-specified value
                f"--sjdbGTFfeatureExon CDS "
                f"--runThreadN 4"
            )
            run_command(command)

        if 'align' in args.pipeline_steps:
            if args.type == 'paired':
                paired_files = get_paired_files(args.input_dir, args.extension)
                for r1, r2 in paired_files:
                    logging.info(f"Step 5: Aligning paired-end data: {os.path.basename(r1)} and {os.path.basename(r2)}")
                    command = f"python3 src/polyAMap.py -i {args.input_dir} -g {indexed_genome_path} -o {args.outputdir} -t paired -e {args.extension} {r1} {r2}"
                    run_command(command)
            elif args.type == 'single':
                fastq_files = get_fastq_files(args.input_dir)
                for fastq_file in fastq_files:
                    logging.info(f"Step 5: Aligning single-end data: {os.path.basename(fastq_file)}")
                    command = f"python3 src/polyAMap.py -i {args.input_dir} -g {indexed_genome_path} -o {args.outputdir} -t single {fastq_file}"
                    run_command(command)

        logging.info("Pipeline execution completed successfully.")
        print("Pipeline execution completed successfully.")

    except Exception as e:
        logging.error(f"Error in main script: {e}")
        sys.exit(1)
