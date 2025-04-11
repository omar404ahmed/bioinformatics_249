#!/usr/bin/env python3
"""
Kraken2 Implementation Script
Implementation for Task 3: Building and using Kraken2 database for classification
"""

import os
import sys
import subprocess
import logging
import re
import time
import platform
from pathlib import Path
import psutil

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("kraken2_task.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("kraken2_task")

class Kraken2Runner:
    """Class to handle Kraken2 database building and classification"""
    
    def __init__(self):
        """Initialize paths and check requirements"""
        # Find kraken2 and kraken2-build paths
        self.kraken2_bin = self._find_executable("kraken2")
        self.kraken_build = self._find_executable("kraken2-build")
        
        # Set directories
        self.db_dir = Path("./kraken2_custom_db")
        self.output_dir = Path("./kraken2_results")
        
        # Create output directories
        self.db_dir.mkdir(exist_ok=True, parents=True)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Performance metrics
        self.start_time = None
        self.end_time = None
        self.peak_memory = 0.0
        
        logger.info("================ KRAKEN2 COMPARISON TASK ================")
        logger.info(f"Using Kraken2 at: {self.kraken2_bin}")
        logger.info(f"Running on: {platform.system()} {platform.release()}")
    
    def _find_executable(self, name):
        """Find path to executable"""
        try:
            path = subprocess.check_output(["which", name], text=True).strip()
            if not path:
                raise FileNotFoundError(f"{name} not found")
            return path
        except (subprocess.SubprocessError, FileNotFoundError):
            logger.error(f"Could not find {name}. Make sure it's installed and in your PATH.")
            sys.exit(1)
    
    def _run_command(self, cmd, desc, log_file=None):
        """Run a command with timing and error handling"""
        logger.info(f"Running: {desc}")
        
        start_time = time.time()
        start_memory = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024  # in MB
        peak_memory = start_memory
        
        # Prepare logging
        if log_file:
            log_path = self.output_dir / log_file
            with open(log_path, "w") as log_file:
                try:
                    process = subprocess.run(
                        cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True
                    )
                    log_file.write(process.stdout)
                    log_file.write(process.stderr)
                    
                    # Record output to log
                    logger.info(f"Command output saved to {log_path}")
                    
                    end_time = time.time()
                    duration = end_time - start_time
                    
                    # Measure memory usage
                    peak_memory = max(peak_memory, psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024)
                    self.peak_memory = max(self.peak_memory, peak_memory)
                    
                    # Record stats to log
                    log_file.write(f"\n\n--- Performance Stats ---\n")
                    log_file.write(f"Execution time: {duration:.2f} seconds\n")
                    log_file.write(f"Peak memory usage: {peak_memory:.2f} MB\n")
                    
                    logger.info(f"Completed: {desc} in {duration:.2f} seconds")
                    return process
                except subprocess.CalledProcessError as e:
                    log_file.write(f"ERROR: {e}\n")
                    log_file.write(e.stdout)
                    log_file.write(e.stderr)
                    logger.error(f"Error executing command: {' '.join(cmd)}")
                    logger.error(e.stderr)
                    return None
        else:
            try:
                process = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
                
                end_time = time.time()
                duration = end_time - start_time
                
                # Measure memory usage
                peak_memory = max(peak_memory, psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024)
                self.peak_memory = max(self.peak_memory, peak_memory)
                
                logger.info(f"Completed: {desc} in {duration:.2f} seconds")
                return process
            except subprocess.CalledProcessError as e:
                logger.error(f"Error executing command: {' '.join(cmd)}")
                logger.error(e.stderr)
                return None
    
    def download_taxonomy(self):
        """Download NCBI taxonomy"""
        logger.info("Downloading NCBI taxonomy...")
        cmd = [
            self.kraken_build,
            "--download-taxonomy",
            "--db", str(self.db_dir)
        ]
        
        return self._run_command(cmd, "Downloading taxonomy", "taxonomy_download.log")
    
    def add_genomes_to_library(self):
        """Add reference genomes to the library"""
        logger.info("Adding genomes to library...")
        
        # List of reference genomes
        genomes = [
            "./reference_genomes/e_coli.fna",
            "./reference_genomes/b_subtilis.fna",
            "./reference_genomes/m_tuberculosis.fna",
            "./reference_genomes/p_aeruginosa.fna",
            "./reference_genomes/s_aureus.fna"
        ]
        
        # Add each genome
        for genome in genomes:
            cmd = [
                self.kraken_build,
                "--add-to-library", genome,
                "--db", str(self.db_dir)
            ]
            
            result = self._run_command(cmd, f"Adding {genome} to library", "add_library.log")
            if not result:
                logger.warning(f"Failed to add {genome}, continuing with remaining genomes")
        
        return True
    
    def build_database(self):
        """Build the Kraken2 database"""
        logger.info("Building Kraken2 database with stricter matching parameters...")
        cmd = [
            self.kraken_build,
            "--build",
            "--threads", "1",
            "--kmer-len", "35",
            "--minimizer-len", "31",
            "--minimizer-spaces", "0",
            "--db", str(self.db_dir)
        ]
        
        return self._run_command(cmd, "Building database", "db_build.log")
    
    def classify_reads(self, input_file, output_prefix, error_type):
        """Classify reads with Kraken2"""
        logger.info(f"Running Kraken2 on {error_type} reads: {input_file}")
        
        output_file = self.output_dir / f"{output_prefix}_results.txt"
        report_file = self.output_dir / f"{output_prefix}_report.txt"
        log_file = f"{output_prefix}_classification.log"
        
        cmd = [
            self.kraken2_bin,
            "--db", str(self.db_dir),
            "--threads", "1",
            "--output", str(output_file),
            "--report", str(report_file),
            "--use-names",
            "--confidence", "1.0",
            "--quick",
            "--no-unclassified-out",
            input_file
        ]
        
        return self._run_command(cmd, f"Classifying {error_type} reads", log_file)
    
    def combine_results(self):
        """Combine results from R1 and R2 files"""
        logger.info("Combining results from all files...")
        
        # Error-free reads
        error_free_r1 = self.output_dir / "error_free_R1_results.txt"
        error_free_r2 = self.output_dir / "error_free_R2_results.txt"
        error_free_combined = self.output_dir / "error_free_combined_results.txt"
        
        with open(error_free_combined, "w") as outfile:
            for file in [error_free_r1, error_free_r2]:
                if file.exists():
                    with open(file, "r") as infile:
                        outfile.write(infile.read())
        
        # Error-containing reads
        error_reads_r1 = self.output_dir / "error_reads_R1_results.txt"
        error_reads_r2 = self.output_dir / "error_reads_R2_results.txt"
        error_reads_combined = self.output_dir / "error_reads_combined_results.txt"
        
        with open(error_reads_combined, "w") as outfile:
            for file in [error_reads_r1, error_reads_r2]:
                if file.exists():
                    with open(file, "r") as infile:
                        outfile.write(infile.read())
        
        # Also combine report files
        error_free_r1_report = self.output_dir / "error_free_R1_report.txt"
        error_free_r2_report = self.output_dir / "error_free_R2_report.txt"
        error_free_combined_report = self.output_dir / "error_free_combined_report.txt"
        
        with open(error_free_combined_report, "w") as outfile:
            for file in [error_free_r1_report, error_free_r2_report]:
                if file.exists():
                    with open(file, "r") as infile:
                        outfile.write(infile.read())
        
        error_reads_r1_report = self.output_dir / "error_reads_R1_report.txt"
        error_reads_r2_report = self.output_dir / "error_reads_R2_report.txt"
        error_reads_combined_report = self.output_dir / "error_reads_combined_report.txt"
        
        with open(error_reads_combined_report, "w") as outfile:
            for file in [error_reads_r1_report, error_reads_r2_report]:
                if file.exists():
                    with open(file, "r") as infile:
                        outfile.write(infile.read())
        
        logger.info("Results combined successfully")
        return True
    
    def extract_metrics(self):
        """Extract and display performance metrics"""
        logger.info("================ PERFORMANCE SUMMARY ================")
        
        # Database building metrics
        db_build_log = self.output_dir / "db_build.log"
        if db_build_log.exists():
            with open(db_build_log, "r") as f:
                log_content = f.read()
                
                # Extract execution time
                time_match = re.search(r"Execution time: (\d+\.\d+) seconds", log_content)
                if time_match:
                    db_time = float(time_match.group(1))
                    logger.info(f"Database Building - Execution time: {db_time:.2f} seconds")
                
                # Extract memory usage
                memory_match = re.search(r"Peak memory usage: (\d+\.\d+) MB", log_content)
                if memory_match:
                    db_memory = float(memory_match.group(1))
                    logger.info(f"Database Building - Peak Memory: {db_memory:.2f} MB")
        
        # Classification metrics - error-free
        error_free_logs = [
            self.output_dir / "error_free_R1_classification.log",
            self.output_dir / "error_free_R2_classification.log"
        ]
        
        total_time = 0
        total_memory = 0
        for log_file in error_free_logs:
            if log_file.exists():
                with open(log_file, "r") as f:
                    log_content = f.read()
                    
                    time_match = re.search(r"Execution time: (\d+\.\d+) seconds", log_content)
                    if time_match:
                        total_time += float(time_match.group(1))
                    
                    memory_match = re.search(r"Peak memory usage: (\d+\.\d+) MB", log_content)
                    if memory_match:
                        total_memory += float(memory_match.group(1))
        
        logger.info(f"Error-free Reads Classification - Execution time: {total_time:.2f} seconds")
        logger.info(f"Error-free Reads Classification - Peak Memory: {total_memory:.2f} MB")
        
        # Classification metrics - error-containing
        error_reads_logs = [
            self.output_dir / "error_reads_R1_classification.log",
            self.output_dir / "error_reads_R2_classification.log"
        ]
        
        total_time = 0
        total_memory = 0
        for log_file in error_reads_logs:
            if log_file.exists():
                with open(log_file, "r") as f:
                    log_content = f.read()
                    
                    time_match = re.search(r"Execution time: (\d+\.\d+) seconds", log_content)
                    if time_match:
                        total_time += float(time_match.group(1))
                    
                    memory_match = re.search(r"Peak memory usage: (\d+\.\d+) MB", log_content)
                    if memory_match:
                        total_memory += float(memory_match.group(1))
        
        logger.info(f"Error-containing Reads Classification - Execution time: {total_time:.2f} seconds")
        logger.info(f"Error-containing Reads Classification - Peak Memory: {total_memory:.2f} MB")
        
        return True
    
    def parse_classification_results(self):
        """Parse and display classification results"""
        logger.info("================ CLASSIFICATION SUMMARY ================")
        
        # Error-free results
        error_free_combined = self.output_dir / "error_free_combined_results.txt"
        if error_free_combined.exists():
            with open(error_free_combined, "r") as f:
                content = f.read()
                lines = content.strip().split('\n')
                total_reads = len(lines)
                classified_reads = sum(1 for line in lines if line.startswith("C\t"))
                unclassified_reads = sum(1 for line in lines if line.startswith("U\t"))
                
                classified_percent = (classified_reads / total_reads * 100) if total_reads > 0 else 0
                unclassified_percent = (unclassified_reads / total_reads * 100) if total_reads > 0 else 0
                
                logger.info("Error-free Reads Results (Combined R1+R2):")
                logger.info(f"Total reads processed: {total_reads}")
                logger.info(f"Classified reads: {classified_reads} ({classified_percent:.2f}%)")
                logger.info(f"Unclassified reads: {unclassified_reads} ({unclassified_percent:.2f}%)")
        
        # Error-containing results
        error_reads_combined = self.output_dir / "error_reads_combined_results.txt"
        if error_reads_combined.exists():
            with open(error_reads_combined, "r") as f:
                content = f.read()
                lines = content.strip().split('\n')
                total_reads = len(lines)
                classified_reads = sum(1 for line in lines if line.startswith("C\t"))
                unclassified_reads = sum(1 for line in lines if line.startswith("U\t"))
                
                classified_percent = (classified_reads / total_reads * 100) if total_reads > 0 else 0
                unclassified_percent = (unclassified_reads / total_reads * 100) if total_reads > 0 else 0
                
                logger.info("Error-containing Reads Results (Combined R1+R2):")
                logger.info(f"Total reads processed: {total_reads}")
                logger.info(f"Classified reads: {classified_reads} ({classified_percent:.2f}%)")
                logger.info(f"Unclassified reads: {unclassified_reads} ({unclassified_percent:.2f}%)")
        
        # Species-specific classification
        species_counts = self.parse_species_classification()
        
        # Format results for output
        logger.info("\nClassification results:")
        logger.info("------------------------")
        logger.info("Classification,Error-free(%),Error-containing(%)")
        
        species_list = ["E. coli", "B. subtilis", "P. aeruginosa", "S. aureus", "M. tuberculosis", "Multi-match (LCA to Bacteria/Gammaproteob.)", "Unclassified"]
        
        for species in species_list:
            error_free_pct = species_counts.get('error_free', {}).get(species, 0)
            error_containing_pct = species_counts.get('error_containing', {}).get(species, 0)
            logger.info(f"{species},{error_free_pct:.2f},{error_containing_pct:.2f}")
        
        # Write formatted results to file
        summary_file = self.output_dir / "classification_summary.csv"
        with open(summary_file, "w") as f:
            f.write("Classification,Error-free(%),Error-containing(%)\n")
            for species in species_list:
                error_free_pct = species_counts.get('error_free', {}).get(species, 0)
                error_containing_pct = species_counts.get('error_containing', {}).get(species, 0)
                f.write(f"{species},{error_free_pct:.2f},{error_containing_pct:.2f}\n")
        
        logger.info("\nPerformance summary:")
        # logger.info(f"Total execution time: {self.end_time - self.start_time:.2f} seconds")
        if self.start_time is not None and self.end_time is not None:
            execution_time = self.end_time - self.start_time
            logger.info(f"Total execution time: {execution_time:.2f} seconds")
        else:
            logger.info("Execution time: Unknown (timing data not available)")
            execution_time = 0  # Default value
        logger.info(f"Peak memory usage: {self.peak_memory:.2f} MB")
        
        # Write performance summary to file
        perf_file = self.output_dir / "performance_summary.txt"
        with open(perf_file, "w") as f:
            f.write("Performance summary:\n")
            if self.start_time is not None and self.end_time is not None:
                execution_time = self.end_time - self.start_time
                f.write(f"Total execution time: {execution_time:.2f} seconds\n")
            else:
                f.write("Total execution time: Unknown (timing data not available)\n")
            f.write(f"Peak memory usage: {self.peak_memory:.2f} MB\n")
        
        logger.info("==================================================")
        logger.info(f"See detailed results in the {self.output_dir} directory")
        
        return True
    
    def parse_species_classification(self):
        """Parse species-specific classification from report files"""
        species_counts = {
            'error_free': {},
            'error_containing': {}
        }
        
        # Extract from error-free combined report
        error_free_report = self.output_dir / "error_free_combined_report.txt"
        if error_free_report.exists():
            with open(error_free_report, "r") as f:
                content = f.read()
                lines = content.split('\n')
                
                for line in lines:
                    if not line.strip():
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 6:
                        continue
                    
                    percentage = float(parts[0])
                    name = parts[5].strip()
                    
                    # Check for each species
                    if "Escherichia coli" in name:
                        species_counts['error_free']["E. coli"] = percentage
                    elif "Bacillus subtilis" in name:
                        species_counts['error_free']["B. subtilis"] = percentage
                    elif "Pseudomonas aeruginosa" in name:
                        species_counts['error_free']["P. aeruginosa"] = percentage
                    elif "Staphylococcus aureus" in name:
                        species_counts['error_free']["S. aureus"] = percentage
                    elif "Mycobacterium tuberculosis" in name:
                        species_counts['error_free']["M. tuberculosis"] = percentage
                    elif "unclassified" in name.lower():
                        species_counts['error_free']["Unclassified"] = percentage
        
        # Extract from error-containing combined report
        error_containing_report = self.output_dir / "error_reads_combined_report.txt"
        if error_containing_report.exists():
            with open(error_containing_report, "r") as f:
                content = f.read()
                lines = content.split('\n')
                
                for line in lines:
                    if not line.strip():
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 6:
                        continue
                    
                    percentage = float(parts[0])
                    name = parts[5].strip()
                    
                    # Check for each species
                    if "Escherichia coli" in name:
                        species_counts['error_containing']["E. coli"] = percentage
                    elif "Bacillus subtilis" in name:
                        species_counts['error_containing']["B. subtilis"] = percentage
                    elif "Pseudomonas aeruginosa" in name:
                        species_counts['error_containing']["P. aeruginosa"] = percentage
                    elif "Staphylococcus aureus" in name:
                        species_counts['error_containing']["S. aureus"] = percentage
                    elif "Mycobacterium tuberculosis" in name:
                        species_counts['error_containing']["M. tuberculosis"] = percentage
                    elif "unclassified" in name.lower():
                        species_counts['error_containing']["Unclassified"] = percentage
        
        # Set defaults for any missing values
        for error_type in ['error_free', 'error_containing']:
            for species in ["E. coli", "B. subtilis", "P. aeruginosa", "S. aureus", "M. tuberculosis", "Multi-match (LCA to Bacteria/Gammaproteob.)", "Unclassified"]:
                if species not in species_counts[error_type]:
                    species_counts[error_type][species] = 0.0
        
        return species_counts
    
    def run_pipeline(self):
        """Run the entire pipeline"""
        self.start_time = time.time()
        
        try:
            # 1. Download taxonomy
            if not self.download_taxonomy():
                logger.error("Failed to download taxonomy")
                return False
            
            # 2. Add genomes to library
            if not self.add_genomes_to_library():
                logger.error("Failed to add all genomes to library")
                # Continue anyway
            
            # 3. Build database
            if not self.build_database():
                logger.error("Failed to build database")
                return False
            
            # 4. Classify reads
            # Error-free reads
            if not self.classify_reads("./sequencing_reads/simulated_reads_no_errors_10k_R1.fastq", "error_free_R1", "error-free"):
                logger.error("Failed to classify error-free R1 reads")
                # Continue anyway
            
            if not self.classify_reads("./sequencing_reads/simulated_reads_no_errors_10k_R2.fastq", "error_free_R2", "error-free"):
                logger.error("Failed to classify error-free R2 reads")
                # Continue anyway
            
            # Error-containing reads
            if not self.classify_reads("./sequencing_reads/simulated_reads_miseq_10k_R1.fastq", "error_reads_R1", "error-containing"):
                logger.error("Failed to classify error-containing R1 reads")
                # Continue anyway
            
            if not self.classify_reads("./sequencing_reads/simulated_reads_miseq_10k_R2.fastq", "error_reads_R2", "error-containing"):
                logger.error("Failed to classify error-containing R2 reads")
                # Continue anyway
            
            # 5. Combine results
            if not self.combine_results():
                logger.error("Failed to combine results")
                return False
            
            # 6. Extract metrics
            if not self.extract_metrics():
                logger.error("Failed to extract metrics")
                return False
            
            # 7. Parse classification results
            if not self.parse_classification_results():
                logger.error("Failed to parse classification results")
                return False
            
            self.end_time = time.time()
            logger.info(f"Pipeline completed in {self.end_time - self.start_time:.2f} seconds")
            return True
        
        except Exception as e:
            logger.exception(f"Error in pipeline: {e}")
            self.end_time = time.time()
            return False

if __name__ == "__main__":
    try:
        # Run the pipeline
        kraken_runner = Kraken2Runner()
        success = kraken_runner.run_pipeline()
        
        # Exit with appropriate code
        if success:
            sys.exit(0)
        else:
            sys.exit(1)
    
    except Exception as e:
        logger.exception(f"Unhandled exception: {e}")
        sys.exit(1)
