#!/usr/bin/env python3
"""
Resource-efficient metagenome read classification using concurrent processing
and streaming approach to minimize memory usage.
"""

import os
import sys
import time
import gc
import psutil
import argparse
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from Bio import SeqIO
from collections import defaultdict, Counter
from tqdm import tqdm
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger('read_classifier')

# Define constants
READ_BATCH_SIZE = 1000  # Process reads in small batches
MAX_MEMORY_PERCENT = 80  # Maximum memory usage before scaling back workers

def measure_memory():
    """Measure current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # Convert to MB

def get_available_memory():
    """Get available system memory in MB."""
    return psutil.virtual_memory().available / 1024 / 1024  # Convert to MB

def is_gzipped(filepath):
    """Check if a file is gzipped by examining its first bytes."""
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def load_genome_chunk(genome_path, chunk_size=50_000_000):
    """
    Load a genome in chunks to reduce memory footprint.
    
    Args:
        genome_path (str): Path to genome FASTA file
        chunk_size (int): Size of each chunk to yield
        
    Yields:
        str: Chunks of genome sequence
    """
    try:
        # Check if file exists
        if not os.path.exists(genome_path):
            logger.error(f"Genome file not found: {genome_path}")
            return
        
        # Open file with gzip if necessary
        if is_gzipped(genome_path):
            f = gzip.open(genome_path, 'rt')
        else:
            f = open(genome_path, 'r')
        
        # Process in chunks
        current_chunk = ""
        for record in SeqIO.parse(f, 'fasta'):
            seq = str(record.seq).upper()
            
            # If adding this sequence would exceed chunk size, yield current chunk
            if len(current_chunk) + len(seq) > chunk_size and current_chunk:
                yield current_chunk
                current_chunk = ""
            
            current_chunk += seq
        
        # Yield the final chunk
        if current_chunk:
            yield current_chunk
            
        f.close()
        
    except Exception as e:
        logger.error(f"Error loading genome {genome_path}: {e}")
        yield ""

def load_reads_batch(fastq_path, start_pos, batch_size):
    """
    Load a batch of reads from a FASTQ file starting at a specific position.
    
    Args:
        fastq_path (str): Path to FASTQ file
        start_pos (int): Position to start reading from (0 for beginning)
        batch_size (int): Number of reads to load
        
    Returns:
        tuple: (reads, new_position) where reads is a list of sequences and
               new_position is the file position after reading the batch
    """
    reads = []
    try:
        # Check if file exists
        if not os.path.exists(fastq_path):
            logger.error(f"Read file not found: {fastq_path}")
            return reads, -1
        
        # Open file with gzip if necessary
        if is_gzipped(fastq_path):
            f = gzip.open(fastq_path, 'rt')
        else:
            f = open(fastq_path, 'r')
        
        # Seek to the starting position if not starting from beginning
        if start_pos > 0:
            f.seek(start_pos)
        
        # Read batch of sequences
        count = 0
        for record in SeqIO.parse(f, 'fastq'):
            reads.append(str(record.seq).upper())
            count += 1
            if count >= batch_size:
                break
        
        # Store current position for next batch
        current_pos = f.tell()
        f.close()
        
        return reads, current_pos
        
    except Exception as e:
        logger.error(f"Error loading reads from {fastq_path}: {e}")
        return reads, -1

def exact_match_worker(reads, genome_chunk, organism):
    """
    Worker function for exact matching of reads against a genome chunk.
    
    Args:
        reads (list): List of read sequences
        genome_chunk (str): Chunk of genome sequence to match against
        organism (str): Name of the organism
        
    Returns:
        dict: Dictionary mapping read indices to match results
    """
    results = {}
    
    for i, read in enumerate(reads):
        # For exact matching, simply check if read is substring of genome
        if read in genome_chunk:
            results[i] = organism
    
    return results

def approx_match_worker(reads, genome_chunk, organism, max_mismatches=1):
    """
    Worker function for approximate matching with up to max_mismatches.
    
    Args:
        reads (list): List of read sequences
        genome_chunk (str): Chunk of genome sequence to match against
        organism (str): Name of the organism
        max_mismatches (int): Maximum allowed mismatches
        
    Returns:
        dict: Dictionary mapping read indices to match results
    """
    results = {}
    
    for i, read in enumerate(reads):
        # Simple sliding window approach for approximate matching
        read_len = len(read)
        found = False
        
        # Check all possible positions in genome chunk
        for j in range(len(genome_chunk) - read_len + 1):
            # Extract potential match
            potential_match = genome_chunk[j:j+read_len]
            
            # Count mismatches
            mismatches = sum(1 for a, b in zip(read, potential_match) if a != b)
            
            if mismatches <= max_mismatches:
                found = True
                break
                
        if found:
            results[i] = organism
    
    return results

def classify_reads(reads, genome_paths, max_workers, match_type='exact', max_mismatches=1):
    """
    Classify reads against reference genomes using parallel processing.
    
    Args:
        reads (list): List of read sequences
        genome_paths (dict): Dictionary mapping organism names to genome paths
        max_workers (int): Maximum number of worker processes
        match_type (str): Type of matching ('exact' or 'approx')
        max_mismatches (int): Maximum mismatches for approximate matching
        
    Returns:
        dict: Dictionary mapping read indices to lists of matching organisms
    """
    all_matches = defaultdict(list)
    
    # Select the appropriate matching function
    if match_type == 'exact':
        match_func = exact_match_worker
    else:
        match_func = partial(approx_match_worker, max_mismatches=max_mismatches)
    
    # Process each genome
    for organism, path in genome_paths.items():
        logger.info(f"Processing genome: {organism}")
        
        # Create chunks of genome to process
        chunk_id = 0
        for genome_chunk in load_genome_chunk(path):
            if not genome_chunk:
                continue
                
            chunk_id += 1
            logger.info(f"Processing chunk {chunk_id} of {organism} ({len(genome_chunk)} bp)")
            
            # Adapt worker count based on system memory
            avail_mem = get_available_memory()
            current_workers = max(1, min(max_workers, int(avail_mem / 1000)))  # 1 worker per ~1GB
            logger.info(f"Using {current_workers} workers (available memory: {avail_mem:.1f} MB)")
            
            # Process reads in parallel for this genome chunk
            futures = []
            
            # Process reads in smaller sub-batches for better memory management
            read_batches = [reads[i:i+READ_BATCH_SIZE] for i in range(0, len(reads), READ_BATCH_SIZE)]
            
            with ProcessPoolExecutor(max_workers=current_workers) as executor:
                # Submit jobs
                for batch_idx, read_batch in enumerate(read_batches):
                    futures.append(
                        executor.submit(match_func, read_batch, genome_chunk, organism)
                    )
                
                # Collect results as they complete
                for batch_idx, future in enumerate(tqdm(
                        as_completed(futures), 
                        total=len(futures),
                        desc=f"Matching against {organism} (chunk {chunk_id})")):
                    try:
                        batch_results = future.result()
                        
                        # Adjust indices based on batch position
                        base_idx = batch_idx * READ_BATCH_SIZE
                        for read_idx, match_org in batch_results.items():
                            global_idx = base_idx + read_idx
                            all_matches[global_idx].append(match_org)
                            
                    except Exception as e:
                        logger.error(f"Error in worker: {e}")
                
                # Force garbage collection after each chunk
                gc.collect()
    
    return all_matches

def analyze_results(read_matches, total_reads):
    """
    Analyze classification results.
    
    Args:
        read_matches (dict): Dictionary mapping read indices to lists of matching organisms
        total_reads (int): Total number of reads
        
    Returns:
        tuple: (unique_matches, multiple_matches, no_matches)
    """
    multiple_matches = sum(1 for matches in read_matches.values() if len(matches) > 1)
    no_matches = total_reads - len(read_matches)
    unique_matches = sum(1 for matches in read_matches.values() if len(matches) == 1)
    
    logger.info(f"\nTotal reads: {total_reads}")
    logger.info(f"Reads with unique match: {unique_matches} ({unique_matches/total_reads*100:.2f}%)")
    logger.info(f"Reads with multiple matches: {multiple_matches} ({multiple_matches/total_reads*100:.2f}%)")
    logger.info(f"Reads with no matches: {no_matches} ({no_matches/total_reads*100:.2f}%)")
    
    return unique_matches, multiple_matches, no_matches

def generate_report(organism_matches, unique_matches, multiple_matches, no_matches, total_reads, match_type):
    """Generate classification report."""
    title = "EXACT MATCHING RESULTS" if match_type == 'exact' else "APPROXIMATE MATCHING RESULTS"
    logger.info(f"\n===== {title} =====")
    
    # Report matches per organism
    logger.info("\nMatches per organism:")
    for org, count in sorted(organism_matches.items(), key=lambda x: x[1], reverse=True):
        logger.info(f"{org}: {count} reads ({count/total_reads*100:.2f}%)")
    
    # Summary statistics
    logger.info("\nSummary Statistics:")
    logger.info(f"Reads with unique match: {unique_matches} ({unique_matches/total_reads*100:.2f}%)")
    logger.info(f"Reads with multiple matches: {multiple_matches} ({multiple_matches/total_reads*100:.2f}%)")
    logger.info(f"Reads with no matches: {no_matches} ({no_matches/total_reads*100:.2f}%)")

def main():
    parser = argparse.ArgumentParser(description="Resource-Efficient Metagenome Classification")
    parser.add_argument('--r1', help='Path to R1 FASTQ file', required=True)
    parser.add_argument('--r2', help='Path to R2 FASTQ file', required=True)
    parser.add_argument('--base_path', help='Base path for genome files', required=True)
    parser.add_argument('--output', help='Output directory for results', default='results')
    parser.add_argument('--match_type', choices=['exact', 'approx'], default='exact',
                       help='Type of matching to perform')
    parser.add_argument('--mismatches', type=int, default=1, 
                       help='Maximum mismatches for approximate matching')
    parser.add_argument('--workers', type=int, default=0,
                       help='Maximum number of worker processes (0 = auto)')
    args = parser.parse_args()
    
    # Auto-determine number of workers if not specified
    if args.workers <= 0:
        args.workers = max(1, psutil.cpu_count(logical=False) - 1)
    
    logger.info(f"Maximum workers: {args.workers}")
    
    # Define genome paths with dynamic path construction for all 5 organisms
    genome_paths = {
        'E. coli K-12 MG1655': os.path.join(args.base_path, 'e_coli/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna'),
        'B. subtilis 168': os.path.join(args.base_path, 'b_subtilis/ncbi_dataset/data/GCA_000009045.1/GCA_000009045.1_ASM904v1_genomic.fna'),
        'P. aeruginosa PAO1': os.path.join(args.base_path, 'p_aeruginosa/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna'),
        'S. aureus NCTC 8325': os.path.join(args.base_path, 's_aureus/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna'),
        'M. tuberculosis H37Rv': os.path.join(args.base_path, 'm_tuberculosis/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna')
    }
    
    # Validate genome paths
    for org, path in genome_paths.items():
        if not os.path.exists(path):
            logger.warning(f"Genome path for {org} not found: {path}")
    
    # Measure initial memory and start time
    initial_memory = measure_memory()
    start_time = time.time()
    
    # Stream-process reads to minimize memory usage
    all_reads = []
    
    # Load reads from R1
    logger.info(f"Loading reads from {args.r1}...")
    pos = 0
    while True:
        batch, pos = load_reads_batch(args.r1, pos, READ_BATCH_SIZE * args.workers)
        if not batch or pos < 0:
            break
        all_reads.extend(batch)
        logger.info(f"Loaded {len(all_reads)} reads from R1 so far")
    
    # Load reads from R2
    logger.info(f"Loading reads from {args.r2}...")
    pos = 0
    while True:
        batch, pos = load_reads_batch(args.r2, pos, READ_BATCH_SIZE * args.workers)
        if not batch or pos < 0:
            break
        all_reads.extend(batch)
        logger.info(f"Loaded {len(all_reads)} reads total (R1+R2)")
    
    total_reads = len(all_reads)
    logger.info(f"Processing {total_reads} reads")
    
    # Perform classification
    read_matches = classify_reads(
        all_reads, 
        genome_paths, 
        args.workers, 
        match_type=args.match_type,
        max_mismatches=args.mismatches
    )
    
    # Count matches per organism
    organism_matches = Counter()
    for matches in read_matches.values():
        for org in matches:
            organism_matches[org] += 1
    
    # Analyze and report results
    unique_matches, multiple_matches, no_matches = analyze_results(read_matches, total_reads)
    generate_report(organism_matches, unique_matches, multiple_matches, no_matches, total_reads, args.match_type)
    
    # Performance metrics
    end_time = time.time()
    final_memory = measure_memory()
    total_time = end_time - start_time
    
    logger.info("\n===== PERFORMANCE METRICS =====")
    logger.info(f"Total execution time: {total_time:.2f} seconds")
    logger.info(f"Peak memory usage: {final_memory - initial_memory:.2f} MB")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Save performance metrics and results
    with open(os.path.join(args.output, f"{args.match_type}_matching_performance.txt"), 'w') as f:
        f.write(f"Total execution time: {total_time:.2f} seconds\n")
        f.write(f"Peak memory usage: {final_memory - initial_memory:.2f} MB\n")
        f.write(f"Total reads processed: {total_reads}\n")
        f.write(f"Reads with unique match: {unique_matches} ({unique_matches/total_reads*100:.2f}%)\n")
        f.write(f"Reads with multiple matches: {multiple_matches} ({multiple_matches/total_reads*100:.2f}%)\n")
        f.write(f"Reads with no matches: {no_matches} ({no_matches/total_reads*100:.2f}%)\n\n")
        f.write("Matches per organism:\n")
        for org, count in sorted(organism_matches.items(), key=lambda x: x[1], reverse=True):
            f.write(f"{org}: {count} reads ({count/total_reads*100:.2f}%)\n")

if __name__ == "__main__":
    main()
