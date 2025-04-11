#!/usr/bin/env python3
"""
Script to download and process SRA samples for Task 3.2 of Kraken2 analysis
- Downloads Standard-8 Kraken2 database
- Downloads and processes 10 SRA samples (5 human gut, 5 wastewater)
- Runs Kraken2 classification on these samples
"""

import os
import sys
import subprocess
import logging
import time
import argparse
from pathlib import Path
import concurrent.futures
import shutil
import requests
import tarfile
import glob

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("sra_processing.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("sra_processing")

class SRAProcessor:
    """Class to handle SRA sample downloading and processing"""
    
    def __init__(self, args):
        """Initialize with command line arguments"""
        self.args = args
        self.workdir = Path(args.workdir).resolve()
        self.db_dir = self.workdir / "kraken2_standard_db"
        self.sra_dir = self.workdir / "sra_data"
        self.output_dir = self.workdir / "kraken2_results"
        
        # SRA accessions to download
        self.human_gut_samples = [
            "SRR11412973",
            "SRR11412976",
            "SRR11412979",
            "SRR11412980",
            "SRR11412984"
        ]
        
        self.wastewater_samples = [
            "SRR21907296",
            "SRR21907303",
            "SRR21907307",
            "SRR21907332",
            "SRR21907330"
        ]
        
        self.all_samples = self.human_gut_samples + self.wastewater_samples
        
        # Create directories
        self.db_dir.mkdir(exist_ok=True, parents=True)
        self.sra_dir.mkdir(exist_ok=True, parents=True)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Performance tracking
        self.start_time = None
        self.end_time = None
    
    def check_tools(self):
        """Check if required tools are installed"""
        tools_to_check = ["kraken2", "prefetch", "fasterq-dump", "wget", "tar"]
        missing_tools = []
        
        for tool in tools_to_check:
            try:
                result = subprocess.run(
                    ["which", tool],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False
                )
                if result.returncode != 0:
                    missing_tools.append(tool)
            except Exception:
                missing_tools.append(tool)
        
        if missing_tools:
            logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            logger.error("Please install these tools before proceeding")
            if "prefetch" in missing_tools or "fasterq-dump" in missing_tools:
                logger.error("SRA Toolkit can be installed via Homebrew: brew install sratoolkit")
            if "kraken2" in missing_tools:
                logger.error("Kraken2 can be installed via Homebrew: brew install kraken2")
            return False
        
        logger.info("All required tools are installed")
        return True
    
    def run_command(self, cmd, desc, check=True):
        """Run a command with proper error handling"""
        logger.info(f"Running: {desc}")
        
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=check
            )
            
            if result.returncode != 0:
                logger.error(f"Command failed: {' '.join(cmd)}")
                logger.error(f"Error: {result.stderr}")
                if check:
                    raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
            else:
                logger.info(f"Successfully completed: {desc}")
            
            return result
        except Exception as e:
            logger.exception(f"Error running command: {' '.join(cmd)}")
            raise
    
    def download_standard_db(self):
        """Download Standard-8 pre-built Kraken2 database"""
        logger.info("Starting download of Kraken2 Standard-8 database")
        
        # Define database URL and output path
        db_url = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20220926.tar.gz"
        db_archive = self.db_dir / "k2_standard_08gb_20220926.tar.gz"
        
        # Check if the database already exists
        db_extracted = self.db_dir / "hash.k2d"
        if db_extracted.exists() and not self.args.force:
            logger.info("Kraken2 Standard-8 database already exists. Skipping download.")
            return True
        
        try:
            # Download the database
            if not db_archive.exists() or self.args.force:
                logger.info(f"Downloading database from {db_url}")
                
                # Check if wget is available
                try:
                    self.run_command(
                        ["wget", "-c", db_url, "-O", str(db_archive)],
                        "Downloading Kraken2 database with wget"
                    )
                except:
                    # Try with requests as a backup
                    logger.info("wget failed, trying with Python requests")
                    response = requests.get(db_url, stream=True)
                    total_size = int(response.headers.get('content-length', 0))
                    block_size = 1024 * 1024  # 1 MB
                    
                    with open(db_archive, 'wb') as f:
                        for i, data in enumerate(response.iter_content(block_size)):
                            progress = (i * block_size) / total_size * 100
                            logger.info(f"Download progress: {progress:.1f}%")
                            f.write(data)
            
            # Extract the database
            logger.info("Extracting database archive")
            with tarfile.open(db_archive, "r:gz") as tar:
                tar.extractall(path=self.db_dir)
            
            # Verify extraction
            if not db_extracted.exists():
                logger.error("Database extraction failed. Required files not found.")
                return False
            
            logger.info("Kraken2 Standard-8 database download and extraction completed")
            return True
            
        except Exception as e:
            logger.exception(f"Error downloading or extracting Kraken2 Standard-8 database: {e}")
            return False
    
    def find_sra_file(self, sra_id):
        """Find the SRA file in the nested directory structure"""
        # Check multiple possible locations for the SRA file
        possible_locations = [
            self.sra_dir / sra_id / f"{sra_id}.sra",                   # Direct location
            self.sra_dir / sra_id / sra_id / f"{sra_id}.sra",          # Nested location
            self.sra_dir / f"{sra_id}.sra"                             # Root location
        ]
        
        # Try to find the SRA file
        for location in possible_locations:
            if location.exists():
                logger.info(f"Found SRA file at: {location}")
                return location
        
        # If no SRA file is found, search for it recursively
        for root, dirs, files in os.walk(self.sra_dir):
            for file in files:
                if file == f"{sra_id}.sra":
                    path = Path(root) / file
                    logger.info(f"Found SRA file at: {path}")
                    return path
        
        return None
    
    def extract_fastq_from_sra(self, sra_id):
        """Extract FASTQ files from SRA file"""
        # First find the SRA file
        sra_file = self.find_sra_file(sra_id)
        if not sra_file:
            logger.error(f"SRA file for {sra_id} not found after recursive search")
            
            # Try to extract directly from accession without finding SRA file
            logger.info(f"Attempting direct FASTQ extraction for {sra_id}")
            sample_dir = self.sra_dir / sra_id
            sample_dir.mkdir(exist_ok=True, parents=True)
            
            try:
                self.run_command(
                    ["fasterq-dump", sra_id, "-O", str(sample_dir), "-t", str(sample_dir), "-e", str(self.args.threads)],
                    f"Direct FASTQ extraction for {sra_id}"
                )
                
                # Check if FASTQ files were created
                fastq_files = list(sample_dir.glob("*.fastq"))
                if fastq_files:
                    logger.info(f"Successfully extracted {len(fastq_files)} FASTQ files for {sra_id}")
                    return True
                else:
                    logger.error(f"No FASTQ files were created for {sra_id}")
                    return False
                    
            except subprocess.CalledProcessError:
                logger.error(f"Failed to extract FASTQ files directly for {sra_id}")
                return False
            
        # Extract FASTQ files from the SRA file
        sample_dir = self.sra_dir / sra_id
        sample_dir.mkdir(exist_ok=True, parents=True)
        
        try:
            self.run_command(
                ["fasterq-dump", str(sra_file), "-O", str(sample_dir), "-t", str(sample_dir), "-e", str(self.args.threads)],
                f"Extracting FASTQ files for {sra_id}"
            )
            
            # Check if FASTQ files were created
            fastq_files = list(sample_dir.glob("*.fastq"))
            if fastq_files:
                logger.info(f"Successfully extracted {len(fastq_files)} FASTQ files for {sra_id}")
                return True
            else:
                logger.error(f"No FASTQ files were created for {sra_id}")
                return False
                
        except subprocess.CalledProcessError:
            logger.error(f"Failed to extract FASTQ files for {sra_id}")
            return False
    
    def download_single_sra(self, sra_id):
        """Download and process a single SRA sample"""
        try:
            sample_dir = self.sra_dir / sra_id
            sample_dir.mkdir(exist_ok=True)
            
            # Check if FASTQ files already exist
            fastq_files = list(sample_dir.glob("*.fastq")) + list(sample_dir.glob("*.fastq.gz"))
            if fastq_files and not self.args.force:
                logger.info(f"FASTQ files for {sra_id} already exist. Skipping download.")
                return True
            
            # Download the SRA file
            logger.info(f"Downloading {sra_id}")
            try:
                self.run_command(
                    ["prefetch", sra_id, "-O", str(self.sra_dir)],
                    f"Downloading {sra_id}"
                )
            except subprocess.CalledProcessError:
                logger.error(f"Failed to download {sra_id} with prefetch")
                return False
            
            # Extract FASTQ files
            if not self.extract_fastq_from_sra(sra_id):
                logger.error(f"Failed to extract FASTQ files for {sra_id}")
                return False
            
            # Compress FASTQ files to save space
            for fastq_file in sample_dir.glob("*.fastq"):
                try:
                    self.run_command(
                        ["gzip", str(fastq_file)],
                        f"Compressing {fastq_file.name}"
                    )
                except subprocess.CalledProcessError:
                    logger.warning(f"Failed to compress {fastq_file.name}")
            
            return True
            
        except Exception as e:
            logger.exception(f"Error processing {sra_id}: {e}")
            return False
    
    def download_sra_samples(self):
        """Download all SRA samples using parallel processing"""
        logger.info(f"Starting download of {len(self.all_samples)} SRA samples")
        
        # Create a mapping file for sample types
        with open(self.output_dir / "sample_types.csv", "w") as f:
            f.write("sample_id,type\n")
            for sample_id in self.human_gut_samples:
                f.write(f"{sample_id},human_gut\n")
            for sample_id in self.wastewater_samples:
                f.write(f"{sample_id},wastewater\n")
        
        # Download samples in parallel
        max_workers = min(self.args.threads, len(self.all_samples))
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_sample = {executor.submit(self.download_single_sra, sample_id): sample_id for sample_id in self.all_samples}
            
            successful_downloads = 0
            for future in concurrent.futures.as_completed(future_to_sample):
                sample_id = future_to_sample[future]
                try:
                    success = future.result()
                    if success:
                        successful_downloads += 1
                except Exception as e:
                    logger.exception(f"Exception downloading {sample_id}: {e}")
        
        logger.info(f"Successfully downloaded {successful_downloads}/{len(self.all_samples)} SRA samples")
        return successful_downloads > 0
    
    def find_fastq_files(self, sample_id):
        """Find FASTQ files for a sample, searching in multiple possible locations"""
        # Base directory for the sample
        sample_dir = self.sra_dir / sample_id
        
        # Look for paired-end files with different naming patterns
        paired_patterns = [
            (f"{sample_id}_1.fastq.gz", f"{sample_id}_2.fastq.gz"),
            (f"{sample_id}.1.fastq.gz", f"{sample_id}.2.fastq.gz")
        ]
        
        # Single-end patterns
        single_patterns = [
            f"{sample_id}.fastq.gz",
            f"{sample_id}.1.fastq.gz"  # Sometimes only R1 is available
        ]
        
        # Check for paired-end files
        for pattern_r1, pattern_r2 in paired_patterns:
            r1_files = list(sample_dir.glob(pattern_r1))
            r2_files = list(sample_dir.glob(pattern_r2))
            
            if r1_files and r2_files:
                logger.info(f"Found paired-end files for {sample_id}: {r1_files[0].name}, {r2_files[0].name}")
                return {"paired": True, "r1": r1_files[0], "r2": r2_files[0]}
        
        # Check for single-end files
        for pattern in single_patterns:
            files = list(sample_dir.glob(pattern))
            if files:
                logger.info(f"Found single-end file for {sample_id}: {files[0].name}")
                return {"paired": False, "file": files[0]}
        
        # If nothing found in the direct sample dir, search recursively
        logger.info(f"Searching recursively for FASTQ files for {sample_id}")
        
        # Recursive search for any FASTQ files
        all_fastq = []
        for root, dirs, files in os.walk(sample_dir):
            for file in files:
                if file.endswith('.fastq.gz') or file.endswith('.fastq'):
                    all_fastq.append(Path(root) / file)
        
        # Try to identify paired files
        r1_files = [f for f in all_fastq if '_1.fastq' in str(f) or '.1.fastq' in str(f)]
        r2_files = [f for f in all_fastq if '_2.fastq' in str(f) or '.2.fastq' in str(f)]
        
        if r1_files and r2_files:
            logger.info(f"Found paired-end files for {sample_id} in recursive search")
            return {"paired": True, "r1": r1_files[0], "r2": r2_files[0]}
        elif all_fastq:
            logger.info(f"Found at least one FASTQ file for {sample_id} in recursive search")
            return {"paired": False, "file": all_fastq[0]}
        
        # No files found
        logger.error(f"No FASTQ files found for {sample_id}")
        return None
    
    def classify_single_sample(self, sample_id):
        """Run Kraken2 classification on a single sample"""
        try:
            # Check if Kraken2 report already exists
            report_file = self.output_dir / f"{sample_id}_report.txt"
            if report_file.exists() and not self.args.force:
                logger.info(f"Kraken2 report for {sample_id} already exists. Skipping classification.")
                return True
            
            # Find FASTQ files
            fastq_info = self.find_fastq_files(sample_id)
            
            if not fastq_info:
                logger.error(f"No FASTQ files found for {sample_id}")
                return False
            
            # Create output files
            output_file = self.output_dir / f"{sample_id}_output.txt"
            
            # Run Kraken2 based on paired or single-end
            if fastq_info["paired"]:
                # Paired-end
                logger.info(f"Running Kraken2 on paired-end reads for {sample_id}")
                
                cmd = [
                    "kraken2",
                    "--db", str(self.db_dir),
                    "--paired",
                    "--output", str(output_file),
                    "--report", str(report_file),
                    "--threads", str(self.args.threads),
                    "--report-minimizer-data",
                    str(fastq_info["r1"]),
                    str(fastq_info["r2"])
                ]
                
            else:
                # Single-end
                logger.info(f"Running Kraken2 on single-end reads for {sample_id}")
                
                cmd = [
                    "kraken2",
                    "--db", str(self.db_dir),
                    "--output", str(output_file),
                    "--report", str(report_file),
                    "--threads", str(self.args.threads),
                    "--report-minimizer-data",
                    str(fastq_info["file"])
                ]
            
            # Run Kraken2
            self.run_command(cmd, f"Classifying {sample_id}")
            
            # Create a clean report file (without comments)
            clean_report = self.output_dir / f"{sample_id}_clean_report.txt"
            with open(report_file, "r") as infile, open(clean_report, "w") as outfile:
                for line in infile:
                    if not line.startswith("#"):
                        outfile.write(line)
            
            return True
            
        except Exception as e:
            logger.exception(f"Error classifying {sample_id}: {e}")
            return False
    
    def classify_sra_samples(self):
        """Run Kraken2 classification on all SRA samples"""
        logger.info(f"Starting classification of {len(self.all_samples)} SRA samples")
        
        # Verify the database exists
        hash_file = self.db_dir / "hash.k2d"
        if not hash_file.exists():
            db_files = list(self.db_dir.glob("*"))
            if len(db_files) > 0:
                # Try to find the database files in a subdirectory
                subdirs = [d for d in self.db_dir.iterdir() if d.is_dir()]
                for subdir in subdirs:
                    if (subdir / "hash.k2d").exists():
                        logger.info(f"Found database in subdirectory: {subdir}")
                        self.db_dir = subdir
                        break
            
            # Check again if the database was found
            if not (self.db_dir / "hash.k2d").exists():
                logger.error("Kraken2 database not found. Please check the database directory.")
                return False
        
        # Process samples sequentially (for better debugging)
        successful_classifications = 0
        for sample_id in self.all_samples:
            try:
                success = self.classify_single_sample(sample_id)
                if success:
                    successful_classifications += 1
            except Exception as e:
                logger.exception(f"Exception classifying {sample_id}: {e}")
        
        logger.info(f"Successfully classified {successful_classifications}/{len(self.all_samples)} SRA samples")
        return successful_classifications > 0
    
    def run_pipeline(self):
        """Run the complete pipeline to download and process SRA samples"""
        self.start_time = time.time()
        
        try:
            # 1. Check required tools
            if not self.check_tools():
                logger.error("Missing required tools. Cannot proceed.")
                return False
            
            # 2. Download Standard-8 Kraken2 database
            if not self.download_standard_db():
                logger.error("Failed to download Kraken2 Standard-8 database")
                return False
            
            # 3. Download SRA samples
            if not self.download_sra_samples():
                logger.error("Failed to download SRA samples")
                return False
            
            # 4. Run Kraken2 classification
            if not self.classify_sra_samples():
                logger.error("Failed to classify SRA samples")
                return False
            
            self.end_time = time.time()
            elapsed_time = self.end_time - self.start_time
            logger.info(f"Pipeline completed in {elapsed_time:.2f} seconds")
            return True
            
        except Exception as e:
            logger.exception(f"Error in pipeline: {e}")
            self.end_time = time.time() if self.start_time else None
            return False

def main():
    """Main function to parse arguments and run the pipeline"""
    parser = argparse.ArgumentParser(description="Download and process SRA samples for Kraken2 analysis")
    
    parser.add_argument("--workdir", type=str, default=".", help="Working directory (default: current directory)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("--force", action="store_true", help="Force redownload even if files exist")
    parser.add_argument("--step", type=str, choices=["download", "extract", "classify", "all"], default="all",
                      help="Run only specific step of the pipeline (default: all)")
    
    args = parser.parse_args()
    
    processor = SRAProcessor(args)
    
    if args.step == "download":
        success = processor.download_sra_samples()
    elif args.step == "extract":
        success = all(processor.extract_fastq_from_sra(sample_id) for sample_id in processor.all_samples)
    elif args.step == "classify":
        success = processor.classify_sra_samples()
    else:  # "all"
        success = processor.run_pipeline()
    
    if success:
        logger.info("Selected tasks completed successfully!")
        return 0
    else:
        logger.error("Failed to complete the selected tasks")
        return 1

if __name__ == "__main__":
    sys.exit(main())
