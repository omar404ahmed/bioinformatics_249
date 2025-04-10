#!/usr/bin/env python3

"""
Simplified BLASTN Analysis Script - Treats each read independently
"""

import os
import subprocess
import time
import shutil
import sys
from collections import defaultdict, Counter
import psutil
import matplotlib.pyplot as plt
import numpy as np
import tempfile

# Define base directory and file paths
BASE_DIR = "/Users/omar/KAUST/Blastn/try1"
GENOMES_DIR = os.path.join(BASE_DIR, "genomes")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
READS_DIR = os.path.join(BASE_DIR, "reads")

# Create results directory if it doesn't exist
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "plots"), exist_ok=True)

# Define genome paths
GENOMES = {
    "E. coli": os.path.join(GENOMES_DIR, "e_coli/ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna"),
    "B. subtilis": os.path.join(GENOMES_DIR, "b_subtilis/ncbi_dataset/data/GCA_000009045.1/GCA_000009045.1_ASM904v1_genomic.fna"),
    "P. aeruginosa": os.path.join(GENOMES_DIR, "p_aeruginosa/ncbi_dataset/data/GCA_000006765.1/GCA_000006765.1_ASM676v1_genomic.fna"),
    "S. aureus": os.path.join(GENOMES_DIR, "s_aureus/ncbi_dataset/data/GCA_000013425.1/GCA_000013425.1_ASM1342v1_genomic.fna"),
    "M. tuberculosis": os.path.join(GENOMES_DIR, "m_tuberculosis/ncbi_dataset/data/GCA_000195955.2/GCA_000195955.2_ASM19595v2_genomic.fna")
}

# Define read file paths - using BOTH R1 and R2
READS = {
    "error_free": [
        os.path.join(READS_DIR, "simulated_reads_no_errors_10k_R1.fasta"),
        os.path.join(READS_DIR, "simulated_reads_no_errors_10k_R2.fasta")
    ],
    "error_containing": [
        os.path.join(READS_DIR, "simulated_reads_miseq_10k_R1.fasta"),
        os.path.join(READS_DIR, "simulated_reads_miseq_10k_R2.fasta")
    ]
}

# Performance tracking variables
start_time = time.time()
memory_tracker = []

def track_memory():
    """Record current memory usage in MB"""
    process = psutil.Process(os.getpid())
    memory_mb = process.memory_info().rss / (1024 * 1024)
    memory_tracker.append((time.time() - start_time, memory_mb))
    return memory_mb

def check_file_exists(file_path):
    """Check if a file exists and print a message if it doesn't"""
    if not os.path.exists(file_path):
        print(f"WARNING: File not found: {file_path}")
        return False
    return True

def check_blast_installation():
    """Check if BLAST+ tools are installed and available"""
    try:
        # Check for blastn
        result = subprocess.run(["which", "blastn"], capture_output=True, text=True)
        if result.returncode != 0:
            print("ERROR: blastn not found in PATH")
            return False
            
        # Check for makeblastdb
        result = subprocess.run(["which", "makeblastdb"], capture_output=True, text=True)
        if result.returncode != 0:
            print("ERROR: makeblastdb not found in PATH")
            return False
        
        return True
    except Exception as e:
        print(f"Error checking BLAST installation: {e}")
        return False

def create_blast_db(genome_files, temp_dir):
    """Create BLAST databases for each genome"""
    db_paths = {}
    
    print("\nCreating BLAST databases...")
    for species, fasta_file in genome_files.items():
        if not check_file_exists(fasta_file):
            continue
            
        # Create a safe name for the database
        species_id = species.replace(" ", "_").replace(".", "_")
        db_path = os.path.join(temp_dir, species_id)
        
        print(f"Creating BLAST database for {species}...")
        
        # Track memory and time
        start_db = time.time()
        pre_mem = track_memory()
        
        # Create the BLAST database
        cmd = ["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", db_path]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            db_paths[species] = db_path
            
            # Track performance
            end_db = time.time()
            post_mem = track_memory()
            print(f"  ✓ Database created successfully in {end_db - start_db:.2f} seconds")
            print(f"  ✓ Memory delta: {post_mem - pre_mem:.2f} MB")
            
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Error creating database: {e}")
            print(f"  ✗ stdout: {e.stdout.decode() if e.stdout else 'None'}")
            print(f"  ✗ stderr: {e.stderr.decode() if e.stderr else 'None'}")
    
    return db_paths

def count_fasta_reads(fasta_file):
    """Count the number of sequences in a FASTA file"""
    if not check_file_exists(fasta_file):
        return 0
        
    count = 0
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception as e:
        print(f"Error counting reads in {fasta_file}: {e}")
        return 0
    
    return count

def run_blastn(query_fasta, db_path, output_file, species, error_type):
    """Run BLASTN search with specific parameters"""
    if not check_file_exists(query_fasta):
        return None
        
    print(f"Running BLASTN: {os.path.basename(query_fasta)} against {species}")
    
    # Track memory and time
    start_blast = time.time()
    pre_mem = track_memory()
    
    # Using parameters to match what we see in the results
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-perc_identity", "100",  # Require 100% identity
        "-strand", "plus",        # Search only the forward strand
        "-max_target_seqs", "1",  # Only report the best match for each read
        "-num_threads", "4"
    ]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Track performance
        end_blast = time.time()
        post_mem = track_memory()
        print(f"  ✓ BLASTN completed in {end_blast - start_blast:.2f} seconds")
        print(f"  ✓ Memory delta: {post_mem - pre_mem:.2f} MB")
        
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Error running BLASTN: {e}")
        print(f"  ✗ stdout: {e.stdout.decode() if e.stdout else 'None'}")
        print(f"  ✗ stderr: {e.stderr.decode() if e.stderr else 'None'}")
        return None

def get_matches_by_read(blast_output):
    """Get mapping of read ID to genome from BLAST output"""
    matches = {}
    if not os.path.exists(blast_output):
        return matches
        
    try:
        with open(blast_output, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 12:
                    read_id = fields[0]
                    matches[read_id] = True  # Just mark that this read matched
    except Exception as e:
        print(f"Error processing BLAST output {blast_output}: {e}")
    
    return matches

def process_reads(read_files, db_paths, error_type, temp_dir):
    """Process all reads against all genomes - simplified approach"""
    print(f"\nProcessing {error_type} reads...")
    
    # Count total reads across all files
    total_reads = 0
    for read_file in read_files:
        file_reads = count_fasta_reads(read_file)
        print(f"Found {file_reads} reads in {os.path.basename(read_file)}")
        total_reads += file_reads
    
    print(f"Total reads to process: {total_reads}")
    
    # Initialize classification counters
    matches_by_genome = {species: 0 for species in db_paths.keys()}
    
    # Track which reads have already been classified
    classified_reads = set()
    
    # Process each genome first (order matters for consistent results)
    for species, db_path in db_paths.items():
        species_matches = set()
        
        # Process each read file against this genome
        for read_file in read_files:
            if not check_file_exists(read_file):
                continue
                
            # Create output file name
            safe_species = species.replace(" ", "_").replace(".", "_")
            read_name = os.path.basename(read_file).replace(".fasta", "")
            output_file = os.path.join(temp_dir, f"{error_type}_{read_name}_vs_{safe_species}.out")
            
            # Run BLASTN
            if run_blastn(read_file, db_path, output_file, species, error_type):
                # Get matching reads
                file_matches = get_matches_by_read(output_file)
                
                # Only count reads not already classified
                new_matches = set(file_matches.keys()) - classified_reads
                species_matches.update(new_matches)
                
        # Update counters
        matches_by_genome[species] = len(species_matches)
        
        # Mark these reads as classified
        classified_reads.update(species_matches)
    
    # Count multi-matches - reads that match multiple genomes
    # This step is skipped because we're processing genomes in order and only
    # counting a read once for the first genome it matches
    multi_matches = 0
    
    # Count unclassified reads
    unclassified = total_reads - sum(matches_by_genome.values()) - multi_matches
    
    # Assemble results
    classification = {
        "total_reads": total_reads,
        "matches_by_genome": matches_by_genome,
        "multi_match": multi_matches,
        "unclassified": unclassified
    }
    
    return classification

def calculate_percentages(classification):
    """Calculate percentages for each category"""
    total_reads = classification["total_reads"]
    percentages = {}
    
    # Genome-specific matches
    for species, count in classification["matches_by_genome"].items():
        percentages[species] = (count / total_reads) * 100
    
    # Multi-match - using proper label from reference
    percentages["Multi-match (LCA to Bacteria/Gammaproteob.)"] = (classification["multi_match"] / total_reads) * 100
    
    # Unclassified
    percentages["Unclassified"] = (classification["unclassified"] / total_reads) * 100
    
    return percentages

def save_classification_summary(error_free_percentages, error_containing_percentages):
    """Save classification summary to CSV file"""
    output_file = os.path.join(RESULTS_DIR, "classification_summary.csv")
    
    with open(output_file, 'w') as f:
        f.write("Classification,Error-free(%),Error-containing(%)\n")
        
        # Write genome-specific percentages
        for species in GENOMES.keys():
            error_free = error_free_percentages.get(species, 0)
            error_containing = error_containing_percentages.get(species, 0)
            f.write(f"{species},{error_free:.2f},{error_containing:.2f}\n")
        
        # Write multi-match and unclassified
        multi_match_key = "Multi-match (LCA to Bacteria/Gammaproteob.)"
        f.write(f"{multi_match_key},{error_free_percentages.get(multi_match_key, 0):.2f},{error_containing_percentages.get(multi_match_key, 0):.2f}\n")
        f.write(f"Unclassified,{error_free_percentages.get('Unclassified', 0):.2f},{error_containing_percentages.get('Unclassified', 0):.2f}\n")
    
    print(f"\nClassification summary saved to {output_file}")
    return output_file

def create_visualization(summary_file):
    """Create visualization of the results"""
    import pandas as pd
    
    # Read the summary file
    df = pd.read_csv(summary_file)
    
    # Create visualization directory
    plot_dir = os.path.join(RESULTS_DIR, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # Create bar charts
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Error-free results
    df.plot(x='Classification', y='Error-free(%)', kind='barh', ax=ax1, color='black')
    ax1.set_title('Classification results for Error-free reads')
    ax1.set_xlabel('Percentage (%)')
    ax1.set_xlim(0, 100)
    ax1.invert_yaxis()  # To match the example image
    
    # Error-containing results
    df.plot(x='Classification', y='Error-containing(%)', kind='barh', ax=ax2, color='black')
    ax2.set_title('Classification results for Error-containing reads')
    ax2.set_xlabel('Percentage (%)')
    ax2.set_xlim(0, 100)
    ax2.invert_yaxis()  # To match the example image
    
    plt.tight_layout()
    viz_file = os.path.join(plot_dir, 'classification_results.png')
    plt.savefig(viz_file)
    plt.close()
    
    # Create table visualization
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        loc='center',
        cellLoc='center',
        colColours=['#f2f2f2']*len(df.columns)
    )
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.5)
    
    plt.title('Table 2: Classification Results for BLAST')
    table_file = os.path.join(plot_dir, 'classification_table.png')
    plt.savefig(table_file, bbox_inches='tight')
    plt.close()
    
    # Create memory usage plot
    mem_times, mem_usages = zip(*memory_tracker)
    
    plt.figure(figsize=(10, 6))
    plt.plot(mem_times, mem_usages, 'o-', color='blue')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage During BLASTN Analysis')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    mem_file = os.path.join(plot_dir, 'memory_usage.png')
    plt.savefig(mem_file)
    plt.close()
    
    print(f"Visualizations saved to {plot_dir}")
    return viz_file, table_file, mem_file

def check_total_percentage(percentages):
    """Check if percentages sum to approximately 100%"""
    total = sum(percentages.values())
    if abs(total - 100) > 1:
        print(f"WARNING: Total percentage {total:.2f}% is not close to 100%")
    else:
        print(f"Total percentage check passed: {total:.2f}%")

def main():
    """Main function to run the entire analysis"""
    print("Starting BLASTN analysis...")
    print(f"Initial memory usage: {track_memory():.2f} MB")
    
    # Check if BLAST+ is installed
    if not check_blast_installation():
        print("ERROR: BLAST+ tools are required but not found")
        return
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp()
    print(f"Using temporary directory: {temp_dir}")
    
    try:
        # Create BLAST databases
        db_paths = create_blast_db(GENOMES, temp_dir)
        
        if not db_paths:
            print("ERROR: Failed to create BLAST databases")
            return
        
        # Process error-free reads
        error_free_classification = process_reads(
            READS["error_free"], 
            db_paths, 
            "error_free", 
            temp_dir
        )
        
        # Process error-containing reads
        error_containing_classification = process_reads(
            READS["error_containing"], 
            db_paths, 
            "error_containing", 
            temp_dir
        )
        
        # Calculate percentages
        error_free_percentages = calculate_percentages(error_free_classification)
        error_containing_percentages = calculate_percentages(error_containing_classification)
        
        # Verify total percentages
        print("\nVerifying classification percentages:")
        print("Error-free percentages:")
        check_total_percentage(error_free_percentages)
        print("Error-containing percentages:")
        check_total_percentage(error_containing_percentages)
        
        # Save summary
        summary_file = save_classification_summary(
            error_free_percentages, 
            error_containing_percentages
        )
        
        # Create visualizations
        viz_file, table_file, mem_file = create_visualization(summary_file)
        
        # Print classification results
        print("\nClassification results:")
        print("------------------------")
        print("Classification,Error-free(%),Error-containing(%)")
        
        for species in GENOMES.keys():
            ef_pct = error_free_percentages.get(species, 0)
            ec_pct = error_containing_percentages.get(species, 0)
            print(f"{species},{ef_pct:.2f},{ec_pct:.2f}")
        
        multi_match_key = "Multi-match (LCA to Bacteria/Gammaproteob.)"
        print(f"{multi_match_key},{error_free_percentages.get(multi_match_key, 0):.2f},{error_containing_percentages.get(multi_match_key, 0):.2f}")
        print(f"Unclassified,{error_free_percentages.get('Unclassified', 0):.2f},{error_containing_percentages.get('Unclassified', 0):.2f}")
        
        # Performance summary
        end_time = time.time()
        total_time = end_time - start_time
        peak_memory = max(mem_usage for _, mem_usage in memory_tracker)
        
        print("\nPerformance summary:")
        print(f"Total execution time: {total_time:.2f} seconds")
        print(f"Peak memory usage: {peak_memory:.2f} MB")
        
    finally:
        # Clean up temporary directory
        print(f"\nCleaning up temporary directory: {temp_dir}")
        shutil.rmtree(temp_dir)
    
    print("\nAnalysis completed successfully")

if __name__ == "__main__":
    main()
