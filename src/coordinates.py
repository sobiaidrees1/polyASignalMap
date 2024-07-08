import os
import argparse
import logging

def convert_gff_to_gtf(gff_file, gtf_file):
    try:
        logging.info(f"Converting GFF to GTF for {gff_file}...")
        with open(gff_file, 'r') as gff_in, open(gtf_file, 'w') as gtf_out:
            for line in gff_in:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) != 9:
                    continue

                if fields[2] in ["CDS", "exon", "gene", "transcript"]:
                    attributes = fields[8]
                    attr_dict = {}
                    for attr in attributes.split(';'):
                        if attr.strip():
                            key_value = attr.strip().split('=')
                            if len(key_value) == 2:
                                key, value = key_value
                                attr_dict[key] = value

                    if "product" in attr_dict:
                        gene_id = attr_dict.pop("product").replace(' ', '')
                        
                        # Remove specific words from gene_id
                        remove_words = ["protein", "polyprotein","poly"]
                        for word in remove_words:
                            gene_id = gene_id.replace(word, '')

                        gene_id = gene_id.strip()
                        
                        # If gene_id becomes empty after removing words, skip
                        if not gene_id:
                            continue

                        gtf_attributes = f'transcript_id "{gene_id}"; gene_id "{gene_id}"; gene_name "{gene_id}";'
                        fields[8] = gtf_attributes
                        gtf_out.write('\t'.join(fields) + '\n')

        logging.info(f"GTF file saved as {gtf_file}")

    except Exception as e:
        logging.error(f"Error converting {gff_file} to GTF: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GFF to GTF.")
    parser.add_argument("-g", "--gfffile", required=True, help="Input GFF file")
    parser.add_argument("-o", "--outputfile", required=True, help="Output GTF file")

    args = parser.parse_args()
    
    logging.basicConfig(filename='coordinates.log', level=logging.INFO)
    
    convert_gff_to_gtf(args.gfffile, args.outputfile)
