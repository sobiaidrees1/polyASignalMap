import argparse
import logging
import subprocess
import os
import shutil

def run_command(command):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e.cmd}")
        logging.error(f"Return code: {e.returncode}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Alignment script for STAR.")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing FASTQ files")
    parser.add_argument("-g", "--genome_dir", required=True, help="Directory containing indexed genome files")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for alignment results")
    parser.add_argument("-t", "--type", choices=['single', 'paired'], required=True, help="Type of sequencing data (single-end or paired-end)")
    parser.add_argument("-e", "--extension", choices=['R2', 'R3'], help="Extension for the second file in paired-end data (R2 or R3)")
    parser.add_argument("-a", "--accession", default="alignment_script", help="Accession number for naming log file")

    # Accept any number of additional arguments for FASTQ files
    parser.add_argument("files", nargs="*", help="FASTQ files to be aligned")

    args = parser.parse_args()

    # Check if the output directory exists, create if not
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Set log file name based on accession number or default name
    log_file = f"{args.accession}_alignment_script.log"
    logging.basicConfig(filename=log_file, level=logging.INFO)

    try:
        if args.type == 'paired':
            if len(args.files) != 2:
                raise ValueError("Expected exactly two input files for paired-end data.")
            
            r1_file = args.files[0]
            r2_file = args.files[1]

            # Extract base filename of R1 file
            base_name = os.path.basename(r1_file)
            base_name = os.path.splitext(base_name)[0].replace('_R1', '')

            # Log the files being processed
            logging.info(f"Processing paired-end files: {r1_file} and {r2_file}")

            command = (
                f"STAR "
                f"--genomeDir {args.genome_dir} "
                f"--readFilesIn {r1_file} {r2_file} "
                f"--outFileNamePrefix {os.path.join(args.output_dir, base_name)} "
                f"--outSAMstrandField intronMotif "
                f"--readFilesCommand zcat "
                f"--outSAMmapqUnique 50 "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--quantMode GeneCounts "
                f"--runThreadN 20 "  # Adjusted to 20 for consistency
                f"--limitBAMsortRAM 1000914045 "  # Adjusted to consistent value
                f"--outSAMattrRGline 'ID:{base_name} SM:{base_name} LB:{base_name} PL:Illumina' "
            )
            logging.info("Starting STAR alignment for paired-end data.")
            run_command(command)

        elif args.type == 'single':
            if len(args.files) != 1:
                raise ValueError("Expected exactly one input file for single-end data.")

            fastq_file = args.files[0]

            # Extract base filename of single-end file
            base_name = os.path.basename(fastq_file)
            base_name = os.path.splitext(base_name)[0]

            # Log the file being processed
            logging.info(f"Processing single-end file: {fastq_file}")

            command = (
                f"STAR "
                f"--genomeDir {args.genome_dir} "
                f"--readFilesIn {fastq_file} "
                f"--outFileNamePrefix {os.path.join(args.output_dir, base_name)} "
                f"--outSAMstrandField intronMotif "
                f"--readFilesCommand zcat "
                f"--outSAMmapqUnique 50 "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--quantMode GeneCounts "
                f"--runThreadN 20 "  # Adjusted to 20 for consistency
                f"--limitBAMsortRAM 1000914045 "  # Adjusted to consistent value
                f"--outSAMattrRGline 'ID:{base_name} SM:{base_name} LB:{base_name} PL:Illumina' "
            )
            logging.info("Starting STAR alignment for single-end data.")
            run_command(command)

        # Create the counts directory inside the output directory and move ReadsPerGene files
        counts_dir = os.path.join(args.output_dir, "counts")
        if not os.path.exists(counts_dir):
            os.makedirs(counts_dir)

        for file in os.listdir(args.output_dir):
            if "ReadsPerGene" in file:
                shutil.move(os.path.join(args.output_dir, file), counts_dir)
                logging.info(f"Moved {file} to {counts_dir}")

    except Exception as e:
        logging.error(f"Error in alignment script: {e}")
        raise
