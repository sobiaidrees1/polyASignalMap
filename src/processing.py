import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Process polyA sequences.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF file")
    parser.add_argument("-o", "--outputdir", required=True, help="Output directory")
    parser.add_argument("-a", "--accession", required=True, help="Accession number")
    parser.add_argument("-u", "--upstream_length", type=int, default=20, help="Length of upstream sequence to extract (default is 20)")
    
    # Adjusting polyA_signals argument to be handled based on main script
    parser.add_argument("-s", "--polyA_signals", nargs='+', help="PolyA signals to search for")

    return parser.parse_args()

def parse_gtf(gtf_file):
    cds_coords = []
    with open(gtf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[2] == 'CDS':
                start = int(fields[3]) - 1  # Convert to 0-based indexing
                end = int(fields[4])
                attributes = fields[8]
                cds_coords.append((start, end, attributes))
    return cds_coords

def parse_attributes(attribute_string):
    attr_dict = {}
    attributes = attribute_string.split(';')
    for attr in attributes:
        if attr.strip():
            key_value = attr.strip().split(maxsplit=1)
            if len(key_value) == 2:
                key, value = key_value
                attr_dict[key] = value.strip('"')
    return attr_dict

def find_polyA_and_upstream(seq, polyA_signals, upstream_length):
    upstream_sequences = []
    current_start = 0
    polyA_detected = False
    polyA_positions = []
    detected_signals = set()  # Set to store detected polyA signals
    
    while current_start < len(seq):
        min_index = len(seq)
        signal = None
        for polyA_signal in polyA_signals:
            signal_index = seq[current_start:].find(polyA_signal)
            if signal_index != -1 and signal_index < min_index:
                min_index = signal_index
                signal = polyA_signal
        
        if signal is None:
            break

        polyA_detected = True
        signal_index = min_index + current_start
        start_index = max(0, signal_index - upstream_length)
        end_index = signal_index + len(signal)
        
        upstream_sequence = seq[start_index:end_index]
        upstream_sequences.append(str(upstream_sequence))  # Convert Seq to str here
        polyA_positions.append(signal_index)
        detected_signals.add(polyA_signal)
        current_start = end_index
    
    if not polyA_detected:
        print(f"Polyadenylation signal not detected in the sequence")
    
    return upstream_sequences, polyA_positions, detected_signals

def process_polyA(fasta_file, gtf_file, output_dir, accession, upstream_length, polyA_signals=None):
    genome = SeqIO.read(fasta_file, "fasta")
    genome_id = genome.id
    
    cds_coords = parse_gtf(gtf_file)
    new_gtf_entries = []
    gene_lengths = []  # List to store gene lengths
    combined_sequence = Seq("")
    
    # Print the polyA signals being used
    if polyA_signals:
        print(f"PolyA signals to search for: {', '.join(polyA_signals)}")
    else:
        print("No specific PolyA signals provided. Using default.")
        polyA_signals = ["AATAAA"]  # Default to AATAAA if polyA_signals is None
    
    print(f"Found {len(cds_coords)} CDS regions in GTF file.")
    
    for idx, (start, end, attributes) in enumerate(cds_coords):
        cds_seq = genome.seq[start:end]
        attr_dict = parse_attributes(attributes)
        gene_name = attr_dict.get("gene_name", f"gene_{idx+1}")
        
        print(f"Processing gene: {gene_name}")
        upstream_seqs, polyA_positions, detected_signals = find_polyA_and_upstream(cds_seq, polyA_signals=polyA_signals, upstream_length=upstream_length)
        
        if upstream_seqs:
            combined_upstream_seq = "".join(upstream_seqs)
            combined_sequence += Seq(combined_upstream_seq)
            
            # Calculate original and final lengths
            original_length = end - start
            final_length = len(combined_upstream_seq)
            
            gene_lengths.append({
                "GeneName": gene_name,
                "OriginalLength": original_length,
                "FinalLength": final_length,
                "DetectedPolyASignals": ', '.join(detected_signals)
            })
            
            # Generate GTF entries
            for polyA_pos, upstream_seq in zip(polyA_positions, upstream_seqs):
                # Adjust the end position to include the length of polyA signal
                original_start = start + polyA_pos - upstream_length
                original_end = start + polyA_pos - 1 + len(upstream_seq)

                new_gtf_entries.append(f"{genome_id}\tGenbank\ttranscript\t{original_start+1}\t{original_end+1}\t.\t+\t.\ttranscript_id \"{gene_name}\"; gene_id \"{gene_name}\"; gene_name \"{gene_name}\";")
                new_gtf_entries.append(f"{genome_id}\tGenbank\texon\t{original_start+1}\t{original_end+1}\t.\t+\t.\ttranscript_id \"{gene_name}\"; gene_id \"{gene_name}\"; gene_name \"{gene_name}\";")
                new_gtf_entries.append(f"{genome_id}\tGenbank\tCDS\t{original_start+1}\t{original_end+1}\t.\t+\t0\ttranscript_id \"{gene_name}\"; gene_id \"{gene_name}\"; gene_name \"{gene_name}\";")
    
    # Write combined sequences to FASTA
    output_fasta = os.path.join(output_dir, f"{accession}_processed.fasta")
    record = SeqRecord(combined_sequence, id=genome_id, description="Combined upstream sequences including polyA signals")
    SeqIO.write([record], output_fasta, "fasta")
    print(f"Combined upstream sequences including polyA signals saved to {output_fasta}")
    
    # Write adjusted GTF entries to GTF file
    output_gtf = os.path.join(output_dir, f"{accession}_processed.gtf")
    with open(output_gtf, 'w') as gtf_out:
        for entry in new_gtf_entries:
            gtf_out.write(entry + '\n')
    print(f"Adjusted GTF file saved to {output_gtf}")
    
    # Write gene lengths to CSV file
    csv_output_file = os.path.join(output_dir, f"{accession}_gene_lengths.csv")
    fieldnames = ["GeneName", "OriginalLength", "FinalLength", "DetectedPolyASignals"]
    with open(csv_output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(gene_lengths)
    print(f"Gene lengths comparison saved to {csv_output_file}")

if __name__ == "__main__":
    args = parse_args()
    
    process_polyA(args.fasta, args.gtf, args.outputdir, args.accession, args.upstream_length, polyA_signals=args.polyA_signals)
