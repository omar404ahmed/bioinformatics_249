import os
import time
import argparse
from Bio import SeqIO
from collections import defaultdict

def get_memory_usage():
    """Get the current memory usage in MB."""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        print("Warning: psutil module not available. Memory usage cannot be measured.")
        return 0

def build_kmer_index(reference_genomes, k=31):
    """Build a k-mer index for reference genomes."""
    kmer_index = defaultdict(lambda: {genome: 0 for genome in reference_genomes})
    total_kmers = 0
    
    for genome_name, genome_path in reference_genomes.items():
        print(f"Processing genome: {genome_name}")
        
        # Read the genome file
        sequences = []
        with open(genome_path, 'r') as file:
            for record in SeqIO.parse(file, 'fasta'):
                sequences.append(str(record.seq).upper())
        
        # Concatenate all sequences
        sequence = ''.join(sequences)
        
        # Extract k-mers and update the index
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' not in kmer:  # Skip k-mers with unknown nucleotides
                kmer_index[kmer][genome_name] += 1
                total_kmers += 1
    
    return kmer_index, total_kmers

def classify_reads(reads_files, kmer_index, k=31, match_threshold=0.1):
    """
    Classify reads based on exact k-mer matches across multiple files.
    
    Parameters:
    - reads_files: List of read files to process
    - kmer_index: Dictionary mapping k-mers to their counts in each genome
    - k: k-mer length
    - match_threshold: Minimum percentage of k-mers in a read that must match a genome for classification
    """
    # Track counts per organism, multi-matches, and unclassified
    classification_results = {
        "E. coli K-12 MG1655": 0,
        "B. subtilis 168": 0,
        "P. aeruginosa PAO1": 0,
        "S. aureus NCTC 8325": 0,
        "M. tuberculosis H37Rv": 0,
        "Multi-match": 0,
        "Unclassified": 0
    }
    
    total_reads = 0
    
    for reads_file in reads_files:
        print(f"Processing reads file: {reads_file}")
        is_error_containing = "miseq" in reads_file
        
        # Read the reads file
        reads = []
        with open(reads_file, 'r') as file:
            for record in SeqIO.parse(file, 'fastq'):
                reads.append(str(record.seq).upper())
        
        # Track total reads
        total_reads += len(reads)
        
        # For each read, extract k-mers and check for matches
        for read in reads:
            # Skip reads shorter than k
            if len(read) < k:
                classification_results["Unclassified"] += 1
                continue
            
            # Extract all k-mers from the read
            read_kmers = []
            for i in range(len(read) - k + 1):
                kmer = read[i:i+k]
                if 'N' not in kmer:
                    read_kmers.append(kmer)
            
            total_kmers = len(read_kmers)
            if total_kmers == 0:
                classification_results["Unclassified"] += 1
                continue
            
            # Count matching k-mers for each genome
            genome_kmer_matches = {genome: 0 for genome in kmer_index[next(iter(kmer_index))].keys()}
            
            for kmer in read_kmers:
                if kmer in kmer_index:
                    for genome, count in kmer_index[kmer].items():
                        if count > 0:
                            genome_kmer_matches[genome] += 1
            
            # Calculate match percentages for each genome
            genome_match_percentages = {}
            for genome, match_count in genome_kmer_matches.items():
                genome_match_percentages[genome] = match_count / total_kmers
            
            # Adjust threshold based on whether the read has errors
            effective_threshold = match_threshold
            # if is_error_containing:
            #     effective_threshold = match_threshold * 0.5  # Lower threshold for error-containing reads
            
            # Find genomes that meet the threshold
            matching_genomes = []
            for genome, percentage in genome_match_percentages.items():
                if percentage >= effective_threshold:
                    matching_genomes.append(genome)
            
            # Classify based on matching results
            if len(matching_genomes) == 0:
                classification_results["Unclassified"] += 1
            elif len(matching_genomes) == 1:
                # Single match
                classification_results[matching_genomes[0]] += 1
            else:
                # Multiple matches
                classification_results["Multi-match"] += 1
    
    # Calculate percentages
    classification_percentages = {}
    for category, count in classification_results.items():
        classification_percentages[category] = (count / total_reads) * 100
    
    return classification_results, classification_percentages, total_reads

def compute_minimizers(sequence, k=31, w=10):
    """Compute minimizers for a sequence."""
    minimizers = set()
    
    # Check if sequence is too short
    if len(sequence) < k or len(sequence) < w + k - 1:
        return minimizers
    
    # Process the sequence in windows
    for i in range(len(sequence) - (w + k - 1) + 1):
        window = sequence[i:i+w+k-1]
        
        # Extract all k-mers in this window
        window_kmers = []
        for j in range(w):
            kmer = window[j:j+k]
            if 'N' not in kmer:
                window_kmers.append(kmer)
        
        if window_kmers:
            # Find the lexicographically smallest k-mer
            minimizer = min(window_kmers)
            minimizers.add(minimizer)
    
    return minimizers

def build_minimizer_index(reference_genomes, k=31, w=10):
    """Build a minimizer-based k-mer index for reference genomes."""
    minimizer_index = defaultdict(lambda: {genome: 0 for genome in reference_genomes})
    total_minimizers = 0
    
    for genome_name, genome_path in reference_genomes.items():
        print(f"Processing genome: {genome_name}")
        
        # Read the genome file
        sequences = []
        with open(genome_path, 'r') as file:
            for record in SeqIO.parse(file, 'fasta'):
                sequences.append(str(record.seq).upper())
        
        # Concatenate all sequences
        sequence = ''.join(sequences)
        
        # Compute minimizers for this sequence
        minimizers = compute_minimizers(sequence, k, w)
        
        # Update the index
        for minimizer in minimizers:
            minimizer_index[minimizer][genome_name] += 1
            total_minimizers += 1
    
    return minimizer_index, total_minimizers

def classify_reads_with_minimizers(reads_files, minimizer_index, k=31, w=10, match_threshold=0.1):
    """
    Classify reads based on minimizer matches across multiple files.
    
    Parameters:
    - reads_files: List of read files to process
    - minimizer_index: Dictionary mapping minimizers to their counts in each genome
    - k: k-mer length
    - w: Window size for minimizers
    - match_threshold: Minimum percentage of minimizers in a read that must match a genome for classification
    """
    # Track counts per organism, multi-matches, and unclassified
    classification_results = {
        "E. coli K-12 MG1655": 0,
        "B. subtilis 168": 0,
        "P. aeruginosa PAO1": 0,
        "S. aureus NCTC 8325": 0,
        "M. tuberculosis H37Rv": 0,
        "Multi-match": 0,
        "Unclassified": 0
    }
    
    total_reads = 0
    
    for reads_file in reads_files:
        print(f"Processing reads file: {reads_file}")
        is_error_containing = "miseq" in reads_file
        
        # Read the reads file
        reads = []
        with open(reads_file, 'r') as file:
            for record in SeqIO.parse(file, 'fastq'):
                reads.append(str(record.seq).upper())
        
        # Track total reads
        total_reads += len(reads)
        
        # For each read, extract minimizers and check for matches
        for read in reads:
            # Skip reads that are too short
            if len(read) < k or len(read) < w + k - 1:
                classification_results["Unclassified"] += 1
                continue
            
            # Compute minimizers for this read
            minimizers = compute_minimizers(read, k, w)
            
            total_minimizers = len(minimizers)
            if total_minimizers == 0:
                classification_results["Unclassified"] += 1
                continue
            
            # Count matching minimizers for each genome
            genome_minimizer_matches = {genome: 0 for genome in minimizer_index[next(iter(minimizer_index))].keys()}
            
            for minimizer in minimizers:
                if minimizer in minimizer_index:
                    for genome, count in minimizer_index[minimizer].items():
                        if count > 0:
                            genome_minimizer_matches[genome] += 1
            
            # Calculate match percentages for each genome
            genome_match_percentages = {}
            for genome, match_count in genome_minimizer_matches.items():
                genome_match_percentages[genome] = match_count / total_minimizers
            
            # Adjust threshold based on whether the read has errors
            effective_threshold = match_threshold
            # if is_error_containing:
            #     effective_threshold = match_threshold * 0.5  # Lower threshold for error-containing reads
            
            # Find genomes that meet the threshold
            matching_genomes = []
            for genome, percentage in genome_match_percentages.items():
                if percentage >= effective_threshold:
                    matching_genomes.append(genome)
            
            # Classify based on matching results
            if len(matching_genomes) == 0:
                classification_results["Unclassified"] += 1
            elif len(matching_genomes) == 1:
                # Single match
                classification_results[matching_genomes[0]] += 1
            else:
                # Multiple matches
                classification_results["Multi-match"] += 1
    
    # Calculate percentages
    classification_percentages = {}
    for category, count in classification_results.items():
        classification_percentages[category] = (count / total_reads) * 100
    
    return classification_results, classification_percentages, total_reads

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Metagenomic classification using k-mers and minimizers')
    parser.add_argument('--reads', nargs='+', required=True, help='FASTQ read files')
    parser.add_argument('--k', type=int, default=31, help='k-mer length')
    parser.add_argument('--w', type=int, default=10, help='Window size for minimizers')
    parser.add_argument('--threshold', type=float, default=0.1, help='Match threshold for classification')
    args = parser.parse_args()
    
    # Define reference genomes
    reference_genomes = {
        "E. coli K-12 MG1655": "genomes/e_coli/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
        "B. subtilis 168": "genomes/b_subtilis/ncbi_dataset/data/GCA_000009045.1/GCA_000009045.1_ASM904v1_genomic.fna",
        "P. aeruginosa PAO1": "genomes/p_aeruginosa/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna",
        "S. aureus NCTC 8325": "genomes/s_aureus/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna",
        "M. tuberculosis H37Rv": "genomes/m_tuberculosis/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"
    }
    
    # Use command-line arguments for reads
    reads_files = args.reads
    k = args.k
    w = args.w  # Window size for minimizers
    match_threshold = args.threshold
    
    # Determine if we're using error-free or error-containing reads
    is_error_free = any("no_errors" in file for file in reads_files)
    dataset_type = "Error-free" if is_error_free else "Error-containing"
    
    # Task 1.1: Build full k-mer index
    print("Task 1.1: Building full k-mer index...")
    start_time = time.time()
    start_memory = get_memory_usage()
    
    kmer_index, total_kmers = build_kmer_index(reference_genomes, k)
    
    end_memory = get_memory_usage()
    end_time = time.time()
    
    print(f"\nResults for Task 1.1:")
    print(f"Full k-mer index built in {end_time - start_time:.2f} seconds")
    print(f"Memory used: {end_memory - start_memory:.2f} MB")
    print(f"Total number of distinct k-mers in index: {len(kmer_index)}")
    print(f"Total number of k-mers (with duplicates): {total_kmers}")
    print(f"Theoretical maximum number of k-mers: {4**k}")
    
    # Task 1.2: Classify reads
    print("\nTask 1.2: Classifying reads with full k-mer index...")
    start_time = time.time()
    
    full_results, full_percentages, total_reads = classify_reads(reads_files, kmer_index, k, match_threshold)
    
    end_time = time.time()
    
    print(f"\nClassification completed in {end_time - start_time:.2f} seconds")
    print(f"Total reads: {total_reads}")
    print(f"Classification results for {dataset_type} reads:")
    print(f"{'Classification':<25} {'Percentage (%)':<15}")
    print(f"{'-'*40}")
    
    # Print results in the desired order
    for category in ["E. coli K-12 MG1655", "B. subtilis 168", "P. aeruginosa PAO1", 
                     "S. aureus NCTC 8325", "M. tuberculosis H37Rv", "Multi-match", "Unclassified"]:
        count = full_results[category]
        percentage = full_percentages[category]
        print(f"{category:<25} {percentage:>14.2f}")
    
    # Task 1.3: Minimizers
    print("\nTask 1.3: Building minimizer-based index...")
    start_time = time.time()
    start_memory = get_memory_usage()
    
    minimizer_index, total_minimizers = build_minimizer_index(reference_genomes, k, w)
    
    end_memory = get_memory_usage()
    end_time = time.time()
    
    print(f"\nResults for Task 1.3:")
    print(f"Minimizer index built in {end_time - start_time:.2f} seconds")
    print(f"Memory used: {end_memory - start_memory:.2f} MB")
    print(f"Total number of distinct minimizers in index: {len(minimizer_index)}")
    print(f"Total number of minimizers (with duplicates): {total_minimizers}")
    print(f"Reduction ratio (compared to full k-mer index): {len(minimizer_index) / len(kmer_index):.2f}")
    
    print("\nClassifying reads with minimizer-based index...")
    start_time = time.time()
    
    minimizer_results, minimizer_percentages, min_total_reads = classify_reads_with_minimizers(
        reads_files, minimizer_index, k, w, match_threshold
    )
    
    end_time = time.time()
    
    print(f"\nClassification completed in {end_time - start_time:.2f} seconds")
    print(f"Total reads: {min_total_reads}")
    print(f"Classification results for {dataset_type} reads:")
    print(f"{'Classification':<25} {'Percentage (%)':<15}")
    print(f"{'-'*40}")
    
    # Print results in the desired order
    for category in ["E. coli K-12 MG1655", "B. subtilis 168", "P. aeruginosa PAO1", 
                     "S. aureus NCTC 8325", "M. tuberculosis H37Rv", "Multi-match", "Unclassified"]:
        count = minimizer_results[category]
        percentage = minimizer_percentages[category]
        print(f"{category:<25} {percentage:>14.2f}")
    
    # Compare full k-mer and minimizer-based results
    print("\nComparison between full k-mer index and minimizer-based index:")
    for category in ["E. coli K-12 MG1655", "B. subtilis 168", "P. aeruginosa PAO1", 
                     "S. aureus NCTC 8325", "M. tuberculosis H37Rv", "Multi-match", "Unclassified"]:
        full_pct = full_percentages[category]
        min_pct = minimizer_percentages[category]
        diff = min_pct - full_pct
        print(f"{category:<25} Full: {full_pct:>6.2f}%  Min: {min_pct:>6.2f}%  Diff: {diff:>+6.2f}%")

if __name__ == "__main__":
    main()
