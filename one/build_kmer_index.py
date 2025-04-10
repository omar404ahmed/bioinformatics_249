#!/usr/bin/env python3
"""
Task 2.1: Build a k-mer index for the 5 reference genomes.
This script extracts all overlapping k-mers from each genome and
builds an index tracking occurrence counts across all genomes.
"""

import os
import sys
import time
import psutil
import pickle
from Bio import SeqIO
from collections import defaultdict

# Constants
K = 31  # k-mer length
BASE_PATH = "genomes"  # Base directory containing the genome files

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

def build_kmer_index(genome_paths, k=31):
    """Build k-mer index for all reference genomes."""
    start_time = time.time()
    print(f"Building k-mer index with k={k}...")
    initial_memory = measure_memory()
    
    # Initialize k-mer index
    # Key: k-mer, Value: list of counts for each organism [e_coli, b_subtilis, p_aeruginosa, s_aureus, m_tuberculosis]
    kmer_index = defaultdict(lambda: [0, 0, 0, 0, 0])
    
    # Process each genome
    total_kmers_per_genome = []
    unique_kmers_per_genome = []
    
    for idx, (name, path) in enumerate(genome_paths.items()):
        print(f"Processing {name}...")
        
        # Load genome
        sequence = ""
        try:
            with open(path, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    sequence += str(record.seq).upper()
        except Exception as e:
            print(f"Error loading genome {name}: {e}")
            continue
            
        # Extract and count k-mers
        kmers = extract_kmers(sequence, k)
        total_kmers_per_genome.append(len(kmers))
        
        # Before adding to index, count unique k-mers in this genome
        unique_kmers = set(kmers)
        unique_kmers_per_genome.append(len(unique_kmers))
        
        # Add to index
        for kmer in kmers:
            kmer_index[kmer][idx] += 1
        
        # Print organism-specific stats
        unique_kmers_in_index = sum(1 for counts in kmer_index.values() if counts[idx] > 0)
        print(f"  Genome size: {len(sequence)} bp")
        print(f"  Total k-mers: {len(kmers)}")
        print(f"  Unique k-mers in this genome: {len(unique_kmers)}")
        print(f"  Unique k-mers in index so far: {unique_kmers_in_index}")
    
    # Calculate memory usage and processing time
    end_time = time.time()
    final_memory = measure_memory()
    
    # Print summary statistics
    print("\nK-mer Index Summary:")
    print(f"Total unique k-mers across all genomes: {len(kmer_index)}")
    print(f"Theoretical maximum possible k-mers: {4**k} k-mers")
    
    # Report shared k-mers statistics
    shared_kmers = defaultdict(int)
    for counts in kmer_index.values():
        genomes_with_kmer = sum(1 for c in counts if c > 0)
        shared_kmers[genomes_with_kmer] += 1
    
    print("\nK-mer sharing statistics:")
    for num_genomes, count in sorted(shared_kmers.items()):
        print(f"  K-mers present in {num_genomes} genomes: {count} ({count/len(kmer_index)*100:.2f}%)")
    
    print(f"\nMemory usage: {final_memory - initial_memory:.2f} MB")
    print(f"Processing time: {end_time - start_time:.2f} seconds")
    
    return kmer_index

def save_index(index, filename):
    """Save the index to a file."""
    print(f"Saving index to {filename}...")
    with open(filename, 'wb') as f:
        pickle.dump(dict(index), f)  # Convert defaultdict to dict for saving
    print(f"Index saved, file size: {os.path.getsize(filename) / (1024*1024):.2f} MB")

if __name__ == "__main__":
    # Define paths to the reference genomes
    genome_paths = {
        'E. coli K-12 MG1655': os.path.join(BASE_PATH, 'e_coli/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna'),
        'B. subtilis 168': os.path.join(BASE_PATH, 'b_subtilis/ncbi_dataset/data/GCA_000009045.1/GCA_000009045.1_ASM904v1_genomic.fna'),
        'P. aeruginosa PAO1': os.path.join(BASE_PATH, 'p_aeruginosa/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna'),
        'S. aureus NCTC 8325': os.path.join(BASE_PATH, 's_aureus/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna'),
        'M. tuberculosis H37Rv': os.path.join(BASE_PATH, 'm_tuberculosis/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna')
    }
    
    # Verify paths exist
    for name, path in genome_paths.items():
        if not os.path.exists(path):
            print(f"Warning: File not found for {name}: {path}")
    
    # Build k-mer index
    kmer_index = build_kmer_index(genome_paths, K)
    
    # Save the index for use in subsequent tasks
    save_index(kmer_index, "kmer_index.pkl")
