# polyA Signal Mapping [polyASignalMap]

## Background
Understanding the regulation of gene expression is fundamental in molecular biology. Polyadenylation is a key post-transcriptional modification where a stretch of adenine nucleotides (A) is added to the 3' end of mRNA molecules. These signals are crucial for mRNA stability, transport from the nucleus to the cytoplasm, and initiation of translation. Polyadenylation signals, typically represented by sequences like "AATAAA" or variants thereof, are recognized by specific protein complexes that cleave the nascent mRNA transcript and add a poly(A) tail. The exact positioning and strength of these signals influence the efficiency and regulation of gene expression. Therefore, accurately mapping these signals within genomic sequences is essential for deciphering gene regulatory mechanisms.

## Key Components
This pipeline automates the workflow for downloading genomic data, processing polyadenylation sequences, and performing sequence alignment. It integrates multiple scripts to facilitate a seamless data processing workflow from start to finish.

## Unique Features

- **Comprehensive PolyA Processing**: Unlike standard pipelines, polyAsignalMap uniquely focuses on detecting polyadenylation signals, and mapping to upstream sequences, rather than complete viral gene sequence.The pipeline works by getting the viral sequence from NCBI, then convert it to individual genes based on the coordinates in the gtf file. This per gene information is then used to identify polyA signal within each gene, and upstream sequence is retried including the polyA signal. These upstream sequence from each gene are then used for mapping to bulk RNAseq data.
- **Integrated Workflow**: Seamless integration of data downloading, format conversion, polyA sequence processing, and alignment simplifies complex viral mapping. Different user friendly options are provided, so that user can adjust mapping according to their data type. Moreover, user can run only a particular step of the pipeline, depending on their need.
<br>![polyAsignal](https://github.com/sobiaidrees1/polyASignalMap/assets/74347135/7d2661a5-eb75-46a9-8a73-47c6ce8d733c) </br>
<span style="font-size: 8px;">**Figure 1.** Workflow of polyAsignal pipeline. The pipeline starts by retrieving the genome and annotation files [GFF] from the NCBI database based on the accession number provided by the user. Genome file is then split into individual gene sequences based on the information in the GFF file, and CDS for each gene are retrieved. The pipeline then looks for polyA signals within individual genes and gets the upstream sequence (20 nuceltodies alongwith polyASignal). If a genes has < defined upstream length i.e. 20 by default, then all nucleotides are retaiend. If a gene has one or more than one polyA Signals, then coordinates of their upstream and polyA sequence are retrieved from the genome, and are saved into a new annotation file. This annotation file is then used to map reads to the reference genome. The pipeline supports both single-end and paired-end sequencing data.</span>
## Usage
<br>polyASignalMap pipeline can be used in two ways.</br> 
<br></br>
**Method # 1:** User can either run it directly :
```bash
python3 polyASignalMap.py -a ACCESSION_NUMBER -i /path/to/input_dir -o /path/to/output_dir -t paired
```
**Command-line Arguments**:
```bash
Arguments:[* = Mandatory]

  -h, --help            show this help message and exit
*  -a ACCESSION, --accession ACCESSION
                        NCBI accession number
   -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output directory (default is current directory)
  -m MAX_READ_LENGTH, --max_read_length MAX_READ_LENGTH
                        Maximum read length for STAR indexing (default is 100)
*  -i INPUT_DIR, --input_dir INPUT_DIR
                        Input directory containing FASTQ files
  -u UPSTREAM_LENGTH, --upstream_length UPSTREAM_LENGTH
                        Length of upstream sequence to extract (default is 20)
 * -t {paired,single}, --type {paired,single}
                        Type of sequencing data (paired-end or single-end)
  -p {download,convert,polyA,align} [{download,convert,polyA,align} ...], --pipeline_steps {download,convert,polyA,align} [{download,convert,polyA,align} ...]
                        Pipeline steps to execute (default: all)
  -x {R2,R3}, --extension {R2,R3}
                        Extension for the second file in paired-end data
                        (default: R2)
  -s POLYA_SIGNALS [POLYA_SIGNALS ...], --polyA_signals POLYA_SIGNALS [POLYA_SIGNALS ...]
                        PolyA signal sequence to detect (default: AATAAA)
  -n GENOMESAINDEXNBASES, --genomeSAindexNbases GENOMESAINDEXNBASES
                        Number of bases to be used for STAR genome indexing
                        (default is 6)
```
<br>Please make sure all required packages and tools are correctly installed.</br>
## Requirements

- Python 3.8 or greater
- STAR aligner
- Required Python packages (listed in `requirements.txt`)

## Installation

1. **Clone the repository**:
    ```bash
    git clone https://github.com/sobiabajwa/polyASignalMap.git
    cd polyASignalMap
    ```

2. **Install dependencies**:
    Ensure you have all necessary dependencies installed. You can install the required Python packages using:
    ```bash
    pip install -r requirements.txt
    ```

3. **Install STAR**:
    Follow the instructions on the [STAR GitHub page](https://github.com/alexdobin/STAR) to install STAR aligner.

### Singularity container

**Method # 2:** A singularity container [polyASignalMap.sif](https://www.dropbox.com/scl/fi/86t4v5pxdvccaolehef2l/polyASignalMap.sif?rlkey=sdz8kzxo5i438gbux0568ba5z&st=0vh1xnze&dl=0) is also provided to avoid installation of required packages. 

### Running the pipeline using singualrity container

To run the pipeline, use the following command:

```bash
singualrity run polyASignalMap.sif -h

Example usage:
singualrity run polyASignalMap.sif  -a ACCESSION_NUMBER -i /path/to/input_dir -o /path/to/output_dir -t paired
```
polyASignalMap is free program, made available under a GNU General Public License. This program comes with ABSOLUTELY NO WARRANTY. See distributed LICENSE for details.

## Any questions?
Please submit an issue to this GitHub repository.
