#!/usr/bin/env python3
"""
Task 2.2: Implement Classification based on k-mer index
This script classifies reads by comparing their k-mers against the previously built k-mer index.
"""

import os
import sys
import time
import pickle
import psutil
from Bio import SeqIO
from collections import defaultdict, Counter

# Constants
K = 31  # k-mer length
ORGANISM_NAMES = [
    'E. coli K-12 MG1655',
    'B. subtilis 168',
    'P. aeruginosa PAO1',
    'S. aureus NCTC 8325',
    'M. tuberculosis H37Rv'
]

def measure_memory():
    """Return the memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def extract_kmers(sequence, k):
    """Extract all k-mers from a sequence."""
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        # Skip k-mers with Ns
        if 'N' not in kmer:
            kmers.append(kmer)
    return kmers

def load_index(filename):
    """Load the index from a file."""
    print(f"Loading index from {filename}...")
    with open(filename, 'rb') as f:
        index = pickle.load(f)
    print(f"Index loaded with {len(index)} unique k-mers")
    return defaultdict(lambda: [0, 0, 0, 0, 0], index)

def load_reads(fastq_path):
    """Load reads from a FASTQ file."""
    reads = []
    try:
        with open(fastq_path, 'r') as f:
            for record in SeqIO.parse(f, 'fastq'):
                reads.append(str(record.seq).upper())
        print(f"Loaded {len(reads)} reads from {fastq_path}")
    except Exception as e:
        print(f"Error loading reads from {fastq_path}: {e}")
    return reads

def classify_reads(reads, kmer_index, k=31):
    """Classify reads based on k-mer matching."""
    start_time = time.time()
    initial_memory = measure_memory()
    print(f"Classifying {len(reads)} reads using k-mer index (k={k})...")
    
    # Initialize counters
    organism_matches = [0, 0, 0, 0, 0]
    read_classifications = []
    
    # Process each read
    for i, read in enumerate(reads):
        if i % 1000 == 0 and i > 0:
            print(f"Processed {i} reads...")
        
        # Extract k-mers from read
        read_kmers = extract_kmers(read, k)
        
        # Count matches to each organism
        matches = [0, 0, 0, 0, 0]
        for kmer in read_kmers:
            if kmer in kmer_index:
                # Add the counts for this k-mer to our running totals
                for j in range(5):
                    if kmer_index[kmer][j] > 0:
                        matches[j] += 1
        
        # Classify based on the organism with most k-mer matches
        if sum(matches) > 0:
            best_match = matches.index(max(matches))
            read_classifications.append(best_match)
            organism_matches[best_match] += 1
        else:
            read_classifications.append(None)  # No match
    
    # Calculate memory usage and processing time
    end_time = time.time()
    final_memory = measure_memory()
    
    # Print classification results
    print("\nClassification Results:")
    print("Matching reads per organism:")
    for i, name in enumerate(ORGANISM_NAMES):
        print(f"{name}:\t{organism_matches[i]} reads")
    print(f"Total classified reads: {sum(organism_matches)}")
    print(f"Unclassified reads: {len(reads) - sum(organism_matches)}")
    print(f"Memory usage: {final_memory - initial_memory:.2f} MB")
    print(f"Classification time: {end_time - start_time:.2f} seconds")
    
    return read_classifications, organism_matches

def classify_reads_advanced(reads, kmer_index, k=31):
    """
    Classify reads with a more sophisticated approach for handling multiple matches.
    Uses a weighted voting system based on k-mer uniqueness.
    """
    start_time = time.time()
    initial_memory = measure_memory()
    print(f"Classifying {len(reads)} reads with advanced method (k={k})...")
    
    # Calculate k-mer specificity weights
    # The more genomes a k-mer appears in, the less specific it is
    print("Calculating k-mer specificity weights...")
    kmer_weights = {}
    for kmer, counts in kmer_index.items():
        # Count how many genomes this k-mer appears in
        genomes_with_kmer = sum(1 for c in counts if c > 0)
        if genomes_with_kmer > 0:
            # Weight is inversely proportional to the number of genomes
            kmer_weights[kmer] = 1.0 / genomes_with_kmer
    
    # Initialize counters
    organism_matches = [0, 0, 0, 0, 0]
    ambiguous_reads = 0
    unclassified_reads = 0
    read_classifications = []
    
    # Process each read
    for i, read in enumerate(reads):
        if i % 1000 == 0 and i > 0:
            print(f"Processed {i} reads...")
        
        # Extract k-mers from read
        read_kmers = extract_kmers(read, k)
        
        # Weighted voting for each organism
        votes = [0.0, 0.0, 0.0, 0.0, 0.0]
        for kmer in read_kmers:
            if kmer in kmer_index:
                weight = kmer_weights.get(kmer, 1.0)
                for j in range(5):
                    if kmer_index[kmer][j] > 0:
                        votes[j] += weight
        
        # Determine classification
        if sum(votes) > 0:
            max_vote = max(votes)
            best_matches = [i for i, v in enumerate(votes) if v > max_vote * 0.9]  # Consider matches within 90% of max
            
            # Check if there's a clear winner or it's ambiguous
            if len(best_matches) == 1:
                best_match = best_matches[0]
                read_classifications.append(best_match)
                organism_matches[best_match] += 1
            else:
                # Ambiguous classification (tied votes)
                read_classifications.append('ambiguous')
                ambiguous_reads += 1
        else:
            # No matches
            read_classifications.append(None)
            unclassified_reads += 1
    
    # Calculate memory usage and processing time
    end_time = time.time()
    final_memory = measure_memory()
    
    # Print classification results
    print("\nAdvanced Classification Results:")
    print("Matching reads per organism:")
    for i, name in enumerate(ORGANISM_NAMES):
        print(f"{name}:\t{organism_matches[i]} reads")
    print(f"Ambiguous classifications: {ambiguous_reads}")
    print(f"Unclassified reads: {unclassified_reads}")
    print(f"Memory usage: {final_memory - initial_memory:.2f} MB")
    print(f"Classification time: {end_time - start_time:.2f} seconds")
    
    return read_classifications, organism_matches

def save_results(classifications, organism_counts, filename):
    """Save classification results to a file."""
    with open(filename, 'w') as f:
        f.write("Classification Results\n")
        f.write("=====================\n\n")
        f.write("Organisms:\n")
        for i, name in enumerate(ORGANISM_NAMES):
            f.write(f"{i}: {name}\n")
        f.write("\nMatching reads per organism:\n")
        for i, name in enumerate(ORGANISM_NAMES):
            f.write(f"{name}:\t{organism_counts[i]} reads\n")
        
        f.write(f"\nTotal classified reads: {sum(count for count in organism_counts if isinstance(count, int))}\n")
        
        # Count ambiguous and unclassified
        ambiguous = classifications.count('ambiguous') if 'ambiguous' in classifications else 0
        unclassified = classifications.count(None)
        
        f.write(f"Ambiguous classifications: {ambiguous}\n")
        f.write(f"Unclassified reads: {unclassified}\n")
    
    print(f"Results saved to {filename}")

if __name__ == "__main__":
    # Check if k-mer index exists
    if not os.path.exists("kmer_index.pkl"):
        print("Error: k-mer index not found. Please run build_kmer_index.py first.")
        sys.exit(1)
    
    # Load k-mer index
    kmer_index = load_index("kmer_index.pkl")
    
    # Define paths to read files
    r1_path = "simulated_reads_no_errors_10k_R1.fastq"
    r2_path = "simulated_reads_no_errors_10k_R2.fastq"
    
    # Verify paths exist
    for path in [r1_path, r2_path]:
        if not os.path.exists(path):
            print(f"Warning: Read file not found: {path}")
    
    # Load reads
    reads_r1 = load_reads(r1_path)
    reads_r2 = load_reads(r2_path)
    reads = reads_r1 + reads_r2
    print(f"Total reads: {len(reads)}")
    
    # Classify reads with basic method
    print("\n--- Basic Classification ---")
    basic_classifications, basic_counts = classify_reads(reads, kmer_index, K)
    save_results(basic_classifications, basic_counts, "basic_classification_results.txt")
    
    # Classify reads with advanced method
    print("\n--- Advanced Classification ---")
    adv_classifications, adv_counts = classify_reads_advanced(reads, kmer_index, K)
    save_results(adv_classifications, adv_counts, "advanced_classification_results.txt")
