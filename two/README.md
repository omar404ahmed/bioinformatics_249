# Genome Assembly and Evaluation

This repository contains implementations of genome assembly algorithms and evaluation scripts for the Bioinformatics Algorithms Assignment 2 (CS249, April 2025). The [report](report2.pdf) outlines the implementation undertaken to achieve the objectives and the strategies employed.

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
  - [De Bruijn Graph Assembler](#de-bruijn-graph-assembler)
  - [Overlap-Layout-Consensus Assembler](#overlap-layout-consensus-assembler)
  - [Example Commands](#example-commands)
- [Evaluation Tools](#evaluation-tools)
  - [QUAST](#quast)
  - [Merqury](#merqury)
  - [Inspector](#inspector)
- [Datasets](#datasets)
- [Expected Output](#expected-output)
- [Troubleshooting](#troubleshooting)
- [Acknowledgements](#acknowledgements)

## Introduction

This project implements two fundamental genome assembly algorithms:

1. **De Bruijn Graph (DBG) Assembly**: Constructs a de Bruijn graph from k-mers and identifies contigs by finding Eulerian paths.
2. **Overlap-Layout-Consensus (OLC) Assembly**: Computes all-vs-all read overlaps, constructs an overlap graph, and generates a consensus sequence for each contig.

The repository also contains scripts for evaluating assemblies using QUAST, Merqury, and Inspector, as well as a Verkko assembly pipeline for a more advanced assembly of the Scincus mitranus (sandfish lizard) genome.

## Installation

### Prerequisites

- Python 3.8+
- Biopython
- NetworkX
- NumPy
- QUAST (for evaluation)
- Bandage (for graph visualization)
- Merqury (for k-mer evaluation and QV score)
- Inspector (for mis-assembly detection)
- Verkko (for advanced assembly)

### Setting up the environment

```bash
# Clone the repository
git clone https://github.com/omar404ahmed/bioinformatics_249/tree/main/two
cd genome-assembly

# Create a conda environment (recommended)
conda create -n assembly python=3.8
conda activate assembly

# Install the required Python packages
pip install biopython networkx numpy

# Install evaluation tools
conda install -c bioconda quast merqury inspector

# Install Bandage for visualization
# Download from https://rrwick.github.io/Bandage/
```

For Ibex cluster users, most tools are pre-installed and can be loaded using:

```bash
module load quast
module load verkko
module load inspector
module load meryl merqury
```

## Directory Structure

```
genome-assembly/
├── README.md
├── debrujin_assembler.py
├── olc_assembler_v2.py
├── ibex_scripts/
│   ├── quast_script.sh
│   ├── merqury_script.sh
│   ├── verkko_assembly.sh
│   ├── verkko1_assembly.sh
│   └── inspec_script.sh
├── data/
│   ├── toy_dataset/
│   │   ├── reference_r.fasta
│   │   ├── reads_r.fastq
│   │   ├── reference_b.fasta
│   │   └── reads_b.fastq
│   ├── mers/
│   │   ├── GCF_000901155.1_ViralProj183710_genomic.fna
│   │   ├── no_error_reads_hiseq_5k.fastq
│   │   ├── no_error_ont_hq_50x.fastq
│   │   ├── reads_hiseq_5k.fastq
└── └── └── ont_hq_50x.fastq
```

## Usage

### De Bruijn Graph Assembler

The De Bruijn Graph assembler takes FASTQ files as input, constructs a de Bruijn graph from k-mers, identifies contigs, and outputs them as a FASTA file.

```bash
python debrujin_assembler.py -k <kmer_size> -i <input_fastq_files> -o <output_fasta> [-g <output_gfa>]
```

Parameters:
- `-k`: k-mer size (required)
- `-i, --input`: Input FASTQ file(s) (required, can specify multiple files)
- `-o, --output`: Output FASTA file path (required)
- `-g, --gfa`: Output GFA file path for Bandage visualization (optional)

### Overlap-Layout-Consensus Assembler

The OLC assembler takes FASTQ files as input, computes all-vs-all read overlaps, constructs an overlap graph, and generates a consensus sequence for each contig.

```bash
python olc_assembler_v2.py --fastq <input_fastq_files> --output <output_prefix> [--min-overlap <min_overlap>] [--min-identity <min_identity>] [--kmer-size <kmer_size>] [--sketch-size <sketch_size>] [--threads <threads>]
```

Parameters:
- `--fastq`: Input FASTQ file(s) (required, can specify multiple files)
- `--output`: Output prefix for files (required)
- `--min-overlap`: Minimum overlap length (default: 20)
- `--min-identity`: Minimum overlap identity (default: 0.9)
- `--kmer-size`: K-mer size for MinHash (default: 15)
- `--sketch-size`: Sketch size for MinHash (default: 100)
- `--threads`: Number of CPU cores to use (default: 1, use 0 for all available cores)

### Example Commands

#### Synthetic Dataset

```bash
# De Bruijn Graph assembly on reads_b.fastq with k=40
python debrujin_assembler.py -k 40 -i data/synthetic/reads_b.fastq -o results/dbg/contigs_b_k40.fasta -g results/dbg/graph_b_k40.gfa

# De Bruijn Graph assembly on reads_r.fastq with k=37
python debrujin_assembler.py -k 35 -i data/synthetic/reads_r.fastq -o results/dbg/contigs_r_k35.fasta -g results/dbg/graph_r_k35.gfa

# De Bruijn Graph assembly on reads_r.fastq with k=47
python debrujin_assembler.py -k 45 -i data/synthetic/reads_r.fastq -o results/dbg/contigs_r_k45.fasta -g results/dbg/graph_r_k45.gfa

# OLC assembly on reads_b.fastq
python olc_assembler_v2.py --fastq data/synthetic/reads_b.fastq --output results/olc/contigs_b --min-overlap 30 --threads 4
```

#### MERS Dataset

```bash
# De Bruijn Graph assembly on error-free HiSeq reads
python debrujin_assembler.py -k 31 -i data/mers/no_error_reads_hiseq_5k.fastq -o results/dbg/mers_hiseq_no_error.fasta

# De Bruijn Graph assembly on HiSeq reads with errors
python ebrujin_assembler.py -k 31 -i data/mers/reads_hiseq_5k.fastq -o results/dbg/mers_hiseq_with_error.fasta

# OLC assembly on error-free ONT reads
python src/olc_assembler_v2.py --fastq data/mers/no_error_ont_hq_50x.fastq --output results/olc/mers_ont_no_error --min-overlap 30 --threads 8

# OLC assembly on ONT reads with errors
python src/olc_assembler_v2.py --fastq data/mers/ont_hq_50x.fastq --output results/olc/mers_ont_with_error --min-overlap 30 --threads 8
```

### Verkko Assembly

For the Scincus mitranus genome assembly:

```bash
# Run Verkko assembly
sbatch scripts/verkko_assembly.sh
```


## Evaluation Tools

### QUAST

QUAST is used to evaluate the quality of genome assemblies.

```bash
# Basic QUAST evaluation
quast.py -o results/evaluation/quast_results -r data/synthetic/reference_r.fasta results/dbg/contigs_r_k35.fasta results/dbg/contigs_r_k45.fasta

# QUAST evaluation for MERS assemblies
quast.py -o results/evaluation/quast_mers -r data/mers/GCF_000901155.1_ViralProj183710_genomic.fna results/dbg/mers_hiseq_no_error.fasta results/olc/mers_ont_no_error.fasta
```

For lizard assembly on Ibex, use the provided script:

```bash
sbatch scripts/quast_script.sh
```

### Merqury

Merqury is used for k-mer based evaluation of assemblies.

```bash
# Run Merqury evaluation
sbatch scripts/merqury_script.sh
```

### Inspector

Inspector is used to detect mis-assemblies.

```bash
# Run Inspector
sbatch scripts/inspec_script.sh
```

## Datasets

All datasets are available at https://bio2vec.cbrc.kaust.edu.sa/data/mowl/cs249_hw2.tar.gz.

### Synthetic Datasets

- Synthetic Dataset 1:
  - `reference_r.fasta`: Reference sequence 1
  - `reads_r.fastq`: First set of simulated reads
  - `reference_b.fasta`: Reference sequence 2
  - `reads_b.fastq`: Second set of simulated reads

- Synthetic Dataset 2 (MERS-CoV):
  - `GCF_000901155.1_ViralProj183710_genomic.fna`: Reference genome
  - Error-free reads: `no_error_reads_hiseq_5k.fastq`, `no_error_ont_hq_50x.fastq`
  - Reads with errors: `reads_hiseq_5k.fastq`, `ont_hq_50x.fastq`

### Real Dataset

Scincus mitranus data from NCBI SRA under study accession SRP563043:
- Hi-C Illumina NovaSeq 6000 paired-end reads (SRR32302809)
- Oxford Nanopore PromethION reads (SRR32302812)
- PacBio Sequel II reads (SRR32302813)
- PacBio Revio reads (SRR32302814)
- RNA-Seq reads for annotation (SRR32302810, SRR32302811)

## Expected Output

### De Bruijn Graph Assembler

```
Loaded 1425 unique k-mers from 1 files
Graph has 1386 nodes and 1425 edges
Found 1 potential contigs
Contig 1: length=1465
============================================================
                    ASSEMBLY METRICS
============================================================
Total assembly length        : 1,465 bp
Number of contigs            : 1
GC content                   : 49.76%
Largest contig               : 1,465 bp
N50                          : 1,465 bp
N90                          : 1,465 bp
L50                          : 1
============================================================
Note: Metrics that require a reference genome are not included.
============================================================

Wrote 1 contigs to results/dbg/contigs_r_k35.fasta
Wrote graph to results/dbg/graph_r_k35.gfa
```

### OLC Assembler

```
Loaded 125 reads from 1 files
Computing sketches for all reads using 4 threads...
Identifying contained reads...
Identified 12 contained reads
Computing overlaps using MinHash filtering with 4 threads...
Built overlap graph with 113 nodes and 245 edges
Finding paths in the overlap graph...
Found 2 paths
Generating layout...
Computing quality-aware consensus sequences...
Contig 1: 29845 bp
Contig 2: 758 bp
Calculating assembly metrics...

ASSEMBLY METRICS:
Total assembly length: 30,603 bp
Number of contigs: 2
GC content: 38.72%
Largest contig: 29,845 bp
N50: 29,845 bp
N90: 758 bp
L50: 1
```

### QUAST Output Example

```
File name: report.txt
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (one contig is here)
Assembly                   contigs_r_k35
# contigs (>= 0 bp)        1             
# contigs (>= 1000 bp)     1             
Total length (>= 0 bp)     1465          
Total length (>= 1000 bp)  1465          
# contigs                  1             
Largest contig             1465          
Total length               1465          
GC (%)                     49.76         
N50                        1465          
N90                        1465          
L50                        1             
L90                        1             
# N's per 100 kbp          0.00          
# mismatches per 100 kbp   0.00          
# indels per 100 kbp       0.00          
Genome fraction (%)        99.93         
Duplication ratio          1.00          
# misassemblies            0             
# misassembled contigs     0             
```

### Merqury Output Example

```
File name: kmer_statistics
Number of 21-mers that are:
  unique             1061469078  (exactly one instance of the kmer is in the input)
  distinct           3038087090  (non-redundant kmer sequences in the input)
  present           30047684002  (...)
  missing         4395008424014  (non-redundant kmer sequences not in the input)

             number of   cumulative   cumulative     presence
              distinct     fraction     fraction   in dataset
frequency        kmers     distinct        total       (1e-6)
--------- ------------ ------------ ------------ ------------
        1   1061469078       0.3494       0.0353     0.000033
        2     86983649       0.3780       0.0411     0.000067
        3     48863965       0.3941       0.0460     0.000100
        4     64592613       0.4154       0.0546     0.000133
        5     87314410       0.4441       0.0691     0.000166
        6    106305236       0.4791       0.0903     0.000200
        7    118430368       0.5181       0.1179     0.000233
        8    122489427       0.5584       0.1506     0.000266
        9    119818168       0.5978       0.1864     0.000300
...
  8240005            1       1.0000       0.9994   274.230952
  8983451            1       1.0000       0.9997   298.973159
  9061269            1       1.0000       1.0000   301.562976
```


```
File Name :merqury_summary
==============================================================
Scincus mitranus (Sandfish) Genome Assembly - Merqury Evaluation
==============================================================
Date: Saturday 10 May 2025 11:31:40 PM +03

K-MER ANALYSIS AND QV SCORE
---------------------------
Assembly file: /ibex/user/ahmedo/backup_assembly/verkko/assembly.fasta
HiFi reads: /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz
K-mer size: 21

QV Score: 42
(QV is a log-scaled measure of the assembly's base-level accuracy)

Estimated error rate: 0.0000631 errors per base
(Lower error rate indicates higher accuracy)
```


### Inspector Output Example

```
Statics of contigs:
Number of contigs       3132
Number of contigs > 10000 bp    3128
Number of contigs >1000000 bp   1155
Total length    3491931429
Total length of contigs > 10000 bp      3491896167
Total length of contigs >1000000bp      2696718955
Longest contig  12917798
Second longest contig length    11796973
N50     2044959
N50 of contigs >1Mbp    2044959

Read to Contig alignment:
Mapping rate /% 99.99
Split-read rate /%      0.74
Depth   8.6146
Mapping rate in large contigs /%        78.14
Split-read rate in large contigs /%     0.59
Depth in large conigs   8.7258

Structural error        47
Expansion       21
Collapse        26
Haplotype switch        0
Inversion       0

Small-scale assembly error /per Mbp     35.0130685888
Total small-scale assembly error        122262
Base substitution       102782
Small-scale expansion   10670
Small-scale collapse    8810

QV      44.1279647878
```


## Troubleshooting

### Common Issues

1. **Memory errors during assembly**: Reduce the dataset size using `seqtk` for testing:
   ```bash
   seqtk sample -s100 input.fastq 0.1 > subsampled.fastq
   ```

2. **Missing dependencies**: Ensure all required packages are installed:
   ```bash
   pip install biopython networkx numpy
   ```

3. **Graph visualization issues**: Make sure Bandage is properly installed:
   ```bash
   # Install Bandage on Ubuntu
   sudo apt-get install bandage
   ```

4. **Cluster job failures**: Check the error logs and adjust resource requirements in the job script:
   ```bash
   # Increase memory for large assemblies
   #SBATCH --mem=256G
   ```

## Acknowledgements

I would like to express my sincere gratitude to:

- Professor Robert Hoehndorf for his wonderful instruction throughout the course and unwavering support
- Teaching Assistants Olga and Fernando for their consistent guidance and helpful feedback
- My classmates Mohammad and Aznaur with whom I had insightful discussions about the implementation of De Bruijn graphs and Verkko assemblies

This project was completed as part of the CS249 Bioinformatics Algorithms course at KAUST in April 2025.
