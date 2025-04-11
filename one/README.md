# CS249 Spring 2025: Metagenomic Classification with k-mer Indexes and Minimizers

This repository contains the implementation for CS249 Spring 2025 Assignment 1, which focuses on metagenomic classification methods using k-mer indexes, minimizers, and comparison with established bioinformatics tools.

## Overview

In metagenomics, DNA sequence reads come from multiple organisms simultaneously, creating a complex mixture of DNA from various sources. This assignment implements different approaches to classify these reads by matching them against reference genomes:

1. **String Matching Approaches**: Implementing efficient algorithms for exact and approximate string matching
2. **K-mer Index**: Building and using k-mer indexes for efficient classification
3. **Minimizers**: Implementing memory-efficient indexing with minimizers
4. **Real-world Applications**: Using Kraken2 for taxonomic classification of both simulated and real metagenomic samples

## Repository Structure

```
.
├── reference_genomes/                  # Reference genome FASTA files
│   ├── e_coli.fna
│   ├── b_subtilis.fna
│   ├── p_aeruginosa.fna
│   ├── s_aureus.fna
│   └── m_tuberculosis.fna
└── sequencing_reads/                   # Simulated reads (both error-free and with errors)
    ├── simulated_reads_no_errors_10k_R1.fastq
    ├── simulated_reads_no_errors_10k_R2.fastq
    ├── simulated_reads_miseq_10k_R1.fastq
    ├── simulated_reads_miseq_10k_R2.fastq
    ├── simulated_reads_no_errors_10k_R1.fasta
    ├── simulated_reads_no_errors_10k_R2.fasta
    ├── simulated_reads_miseq_10k_R1.fasta
    └── simulated_reads_miseq_10k_R2.fasta
```

## Dependencies

To run the scripts in this repository, you'll need the following:

### Basic Dependencies
- Python 3.8+
- Biopython (`pip install biopython`)
- NumPy (`pip install numpy`)
- Pandas (`pip install pandas`)
- Matplotlib (`pip install matplotlib`)
- Seaborn (`pip install seaborn`)
- tqdm (`pip install tqdm`)
- psutil (`pip install psutil`)
- SciPy (`pip install scipy`)
- scikit-learn (`pip install scikit-learn`)

### External Tools
- BLAST+ suite (`blastn`, `makeblastdb`)
- Kraken2 and KrakenTools
- SRA Toolkit (`prefetch`, `fastq-dump`, `fasterq-dump`)

### Installation on Ubuntu/Debian

```bash
# Install Python packages
pip install biopython numpy pandas matplotlib seaborn tqdm psutil scipy scikit-learn

# Install BLAST+
sudo apt-get update
sudo apt-get install ncbi-blast+

# Install Kraken2
sudo apt-get install kraken2

# Install SRA Toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
export PATH=$PATH:$PWD/sratoolkit.*/bin
```

### Installation on macOS

```bash
# Install Python packages
pip install biopython numpy pandas matplotlib seaborn tqdm psutil scipy scikit-learn

# Install BLAST+ and others with Homebrew
brew install blast
brew install kraken2
brew install sratoolkit
```

## Dataset Preparation

Before running the code, you need to download the reference genomes and prepare the dataset:

```bash
# Create directories
mkdir -p reference_genomes sequencing_reads

# Download reference genomes from NCBI
# You can retrieve them based on the accession numbers provided in the assignment:
# E. coli K-12 MG1655 (GCF_000005845.2)
# B. subtilis 168 (GCF_000009045.1)
# P. aeruginosa PAO1 (GCF_000006765.1)
# S. aureus NCTC 8325 (GCF_000013425.1)
# M. tuberculosis H37Rv (GCF_000195955.2)

# Download the simulated reads from the provided repository
wget https://github.com/bio-ontology-research-group/cs249-2025/raw/main/project1-data/simulated_reads_no_errors_10k_R1.fastq -O sequencing_reads/simulated_reads_no_errors_10k_R1.fastq
wget https://github.com/bio-ontology-research-group/cs249-2025/raw/main/project1-data/simulated_reads_no_errors_10k_R2.fastq -O sequencing_reads/simulated_reads_no_errors_10k_R2.fastq
wget https://github.com/bio-ontology-research-group/cs249-2025/raw/main/project1-data/simulated_reads_miseq_10k_R1.fastq -O sequencing_reads/simulated_reads_miseq_10k_R1.fastq
wget https://github.com/bio-ontology-research-group/cs249-2025/raw/main/project1-data/simulated_reads_miseq_10k_R2.fastq -O sequencing_reads/simulated_reads_miseq_10k_R2.fastq
```

## Task 1: Metagenome Classification by String Matching

### Task 1.1: Multiple Matches Strategy

For handling reads that match multiple organisms, a reasonable strategy is to:

1. **Track all matches**: Store all organisms that match each read
2. **Apply confidence thresholds**: Require a minimum match length or quality
3. **Hierarchical assignment**: Use the most specific taxonomic level possible
4. **LCA (Lowest Common Ancestor)**: Assign reads to the LCA of all matching organisms
5. **Report ambiguity**: Flag reads with multiple matches separately

### Task 1.2: Exact Matching Implementation

The repository includes two implementations for efficient exact matching:

#### Memory-Efficient FM-Index (`memory_efficient_fm.py`)

This implementation provides a memory-efficient FM-Index structure for string matching:

```bash
cd task1_string_matching
python memory_efficient_fm.py
```

Key features of this implementation:
- Uses the Burrows-Wheeler Transform (BWT) for efficient string matching
- Implements FM-Index with sampled suffix arrays to reduce memory footprint
- Provides functions for exact matching and matching with mismatches

#### Resource-Efficient Matching (`resource_efficient_matching.py`)

This script implements a parallelized approach to match reads against reference genomes:

```bash
cd task1_string_matching
python resource_efficient_matching.py --r1 ../sequencing_reads/simulated_reads_no_errors_10k_R1.fastq \
                                     --r2 ../sequencing_reads/simulated_reads_no_errors_10k_R2.fastq \
                                     --base_path ../reference_genomes \
                                     --match_type exact \
                                     --output results
```

Key features:
- Processes genomes in chunks to minimize memory usage
- Uses parallel processing to accelerate matching
- Dynamically adapts the number of workers based on available system resources
- Handles both exact and approximate matching

### Task 1.4: Comparison with BLASTN

The `blastn_script.py` script provides a comparison with BLAST for read classification:

```bash
cd task1_string_matching
python blastn_script.py
```

This script:
- Creates BLAST databases for each reference genome
- Runs BLASTN with parameters optimized for read matching
- Processes results to classify reads by organism
- Compares performance metrics (time, memory, accuracy) with the custom implementations
- Generates visualizations of the classification results

## Task 2: Metagenomic classification by k-mer index

### Task 2.1: Build the k-mer Index

The script `build_kmer_index.py` constructs a k-mer index from the reference genomes:

```bash
cd task2_kmer_index
python build_kmer_index.py
```

Key features:
- Extracts all overlapping k-mers (k=31 by default) from each genome
- Builds an index recording the occurrence counts of each k-mer across all genomes
- Tracks memory usage and processing time
- Saves the index for subsequent use in classification

Expected output includes:
- A pickle file (`kmer_index.pkl`) containing the k-mer index
- Summary statistics on unique k-mers, memory usage, and processing time
- Analysis of k-mer sharing patterns across genomes

### Task 2.2: Implement Classification

The script `classify_reads.py` uses the k-mer index to classify reads:

```bash
cd task2_kmer_index
python classify_reads.py
```

Key features:
- Loads the pre-built k-mer index
- Extracts k-mers from each read and matches against the index
- Implements both basic and advanced classification strategies
- Advanced strategy uses weighted voting based on k-mer uniqueness
- Handles multiple matches and ambiguous classifications
- Reports classification accuracy and performance metrics

### Task 2.3: Minimizers

The script `kmer-minimizer-implementation.py` implements a minimizer-based approach to reduce memory requirements:

```bash
cd task2_kmer_index
python kmer-minimizer-implementation.py --reads ../sequencing_reads/simulated_reads_no_errors_10k_R1.fastq ../sequencing_reads/simulated_reads_no_errors_10k_R2.fastq --k 31 --w 10 --threshold 0.1
```

Key features:
- Implements the minimizer scheme to select representative k-mers
- Compares full k-mer index with minimizer-based index
- Analyzes memory reduction and impact on classification accuracy
- Evaluates the trade-off between memory efficiency and classification performance

Expected output includes:
- Comparative analysis between full k-mer index and minimizer-based index
- Memory usage statistics and reduction ratio
- Classification results for both approaches

## Task 3: Real-world data and tools

### Task 3.1-3.2: Kraken2 Implementation and Analysis

The repository includes several scripts for working with Kraken2 and real-world metagenomic data:

#### Building and Using a Custom Kraken2 Database

```bash
cd task3_kraken2
python kraken2_implementation.py
```

This script:
- Downloads NCBI taxonomy
- Adds reference genomes to the Kraken2 library
- Builds a custom Kraken2 database
- Classifies the simulated reads (both error-free and with errors)
- Generates classification reports and performance metrics

#### Processing Real-world SRA Samples

For downloading and processing real metagenomic samples from SRA:

```bash
cd task3_kraken2
python download_process_sra.py --workdir ./kraken2_workdir --threads 4
```

Or using the simplified version:

```bash
cd task3_kraken2
python simple_process_sra.py --workdir ./kraken2_workdir --threads 4
```

These scripts:
- Download the Standard-8 pre-built Kraken2 database
- Download and extract FASTQ files from the specified SRA accessions
- Run Kraken2 classification on these samples
- Generate reports for further analysis

#### Analyzing Metagenomic Data

The script `metagenomic_analysis.py` provides tools for analyzing and visualizing the Kraken2 results:

```bash
cd task3_kraken2
python metagenomic_analysis.py --results-dir ./kraken2_workdir/kraken2_results
```

Key features:
- Creates taxonomic profile matrices from Kraken2 reports
- Performs PCA to visualize sample clustering
- Conducts hierarchical clustering to identify sample relatedness
- Generates heatmaps of taxonomic abundance
- Determines if samples can be separated by environmental source

## Expected Results

### Task 1 Results
- Classification results for both error-free and error-containing reads
- Performance metrics (time, memory) for custom implementations and BLAST
- Comparative analysis showing accuracy and efficiency trade-offs

### Task 2 Results
- Complete k-mer index with statistics on unique k-mers across genomes
- Classification reports showing organism matches for each read set
- Analysis of minimizer approach, demonstrating memory reduction
- Comparison of classification accuracy between full index and minimizer-based approach

### Task 3 Results
- Kraken2 classification results for simulated reads
- Analysis of real metagenomic samples from human gut and wastewater environments
- Visualization showing clear separation between environmental sources
- Performance comparison between Kraken2 and custom implementations

## Troubleshooting

### Common Issues

1. **Missing dependencies**
   - Ensure all Python packages are installed with `pip install -r requirements.txt`
   - Verify external tools (BLAST+, Kraken2, SRA Toolkit) are installed and in PATH

2. **Memory errors**
   - For large genomes, try increasing the chunk size in resource_efficient_matching.py
   - Use minimizer approach for memory-constrained environments

3. **SRA download issues**
   - If prefetch/fasterq-dump fails, try the simplified script with fastq-dump
   - Ensure your internet connection is stable for large downloads

4. **Performance optimization**
   - Adjust the number of threads based on your system specifications
   - For faster processing, consider subsetting the data as mentioned in the assignment instructions

### Getting Help

If you encounter issues not covered here, please:
1. Check the error messages for specific guidance
2. Review the script code for configuration options


## References

- BLAST: Altschul SF, et al. (1990). Basic local alignment search tool. J Mol Biol, 215(3):403-10.
- Kraken2: Wood DE, et al. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1):257.
- Minimizers: Roberts M, et al. (2004). Reducing storage requirements for biological sequence comparison. Bioinformatics, 20(18):3363-9.
- FM-Index: Ferragina P, Manzini G. (2000). Opportunistic data structures with applications. Proceedings of the 41st Symposium on Foundations of Computer Science.

## License

This project is released under the MIT License.

