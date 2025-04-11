#!/usr/bin/env python3
"""
Simplified script to process SRA files and run Kraken2
Uses fastq-dump instead of fasterq-dump for better reliability
"""

import os
import sys
import subprocess
import logging
import time
import argparse
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("sra_processing_simple.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("sra_processing_simple")

class SimpleSRAProcessor:
    """Simplified class to handle SRA processing"""
    
    def __init__(self, args):
        """Initialize with command line arguments"""
        self.args = args
        self.workdir = Path(args.workdir).resolve()
        self.sra_dir = self.workdir / "sra_data"
        self.db_dir = self.workdir / "kraken2_standard_db"
        self.output_dir = self.workdir / "kraken2_results"
        
        # SRA accessions
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
        
        # Create output directory
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
    def run_command(self, cmd, desc, check=False):
        """Run a command with proper error handling"""
        logger.info(f"Running: {desc}")
        try:
            # Print the full command for debugging
            cmd_str = " ".join(str(x) for x in cmd)
            logger.info(f"Command: {cmd_str}")
            
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=check
            )
            
            # Log command output
            if result.stdout:
                logger.info(f"Command stdout: {result.stdout[:1000]}...")
            if result.stderr:
                logger.info(f"Command stderr: {result.stderr[:1000]}...")
                
            if result.returncode != 0:
                logger.warning(f"Command returned non-zero exit status: {result.returncode}")
            else:
                logger.info(f"Successfully completed: {desc}")
                
            return result
        except Exception as e:
            logger.exception(f"Error running command: {e}")
            return None
            
    def process_sra_file(self, accession):
        """Process a single SRA accession using fastq-dump"""
        logger.info(f"Processing SRA accession: {accession}")
        
        # Skip if FASTQ already exists in output directory
        output_file = self.output_dir / f"{accession}_report.txt"
        if output_file.exists() and not self.args.force:
            logger.info(f"Output file {output_file} already exists. Skipping.")
            return True
            
        # Use prefetch and fastq-dump directly on the accession
        try:
            # Run fastq-dump directly on the accession
            fastq_cmd = [
                "fastq-dump", 
                "--split-3",           # Split into R1/R2 for paired-end
                "--outdir", str(self.output_dir),
                accession
            ]
            
            # Run the command
            result = self.run_command(fastq_cmd, f"Extracting FASTQ for {accession}")
            if result and result.returncode == 0:
                logger.info(f"Successfully extracted FASTQ for {accession}")
            else:
                logger.warning(f"FASTQ extraction may have failed for {accession}")
                
            # Check if FASTQ files were created
            r1_file = self.output_dir / f"{accession}_1.fastq"
            r2_file = self.output_dir / f"{accession}_2.fastq"
            
            if r1_file.exists() and r2_file.exists():
                logger.info(f"Found paired-end FASTQ files for {accession}")
                kraken_cmd = [
                    "kraken2",
                    "--db", str(self.db_dir),
                    "--paired",
                    "--output", str(self.output_dir / f"{accession}_output.txt"),
                    "--report", str(self.output_dir / f"{accession}_report.txt"),
                    "--threads", str(self.args.threads),
                    str(r1_file),
                    str(r2_file)
                ]
            elif r1_file.exists():
                logger.info(f"Found single-end FASTQ file for {accession}")
                kraken_cmd = [
                    "kraken2",
                    "--db", str(self.db_dir),
                    "--output", str(self.output_dir / f"{accession}_output.txt"),
                    "--report", str(self.output_dir / f"{accession}_report.txt"),
                    "--threads", str(self.args.threads),
                    str(r1_file)
                ]
            else:
                logger.error(f"No FASTQ files found for {accession}")
                return False
                
            # Run Kraken2
            kraken_result = self.run_command(kraken_cmd, f"Running Kraken2 on {accession}")
            
            # Create a clean report
            if output_file.exists():
                clean_report = self.output_dir / f"{accession}_clean_report.txt"
                with open(output_file, "r") as infile, open(clean_report, "w") as outfile:
                    for line in infile:
                        if not line.startswith("#"):
                            outfile.write(line)
                            
                logger.info(f"Created clean report for {accession}")
                            
            # Clean up FASTQ files to save space
            if not self.args.keep_fastq:
                for fastq_file in [r1_file, r2_file]:
                    if fastq_file.exists():
                        try:
                            fastq_file.unlink()
                            logger.info(f"Removed {fastq_file}")
                        except Exception as e:
                            logger.warning(f"Could not remove {fastq_file}: {e}")
            
            return True
            
        except Exception as e:
            logger.exception(f"Error processing {accession}: {e}")
            return False
            
    def process_all_samples(self):
        """Process all SRA samples"""
        logger.info(f"Processing {len(self.all_samples)} SRA samples")
        
        # Create a mapping file for sample types
        with open(self.output_dir / "sample_types.csv", "w") as f:
            f.write("sample_id,type\n")
            for sample_id in self.human_gut_samples:
                f.write(f"{sample_id},human_gut\n")
            for sample_id in self.wastewater_samples:
                f.write(f"{sample_id},wastewater\n")
        
        successful = 0
        for sample_id in self.all_samples:
            try:
                if self.process_sra_file(sample_id):
                    successful += 1
            except Exception as e:
                logger.exception(f"Error processing {sample_id}: {e}")
        
        logger.info(f"Successfully processed {successful}/{len(self.all_samples)} samples")
        return successful > 0
        
def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Simple script to process SRA files and run Kraken2")
    parser.add_argument("--workdir", type=str, default=".", help="Working directory (default: current directory)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for Kraken2 (default: 4)")
    parser.add_argument("--force", action="store_true", help="Force reprocessing even if output exists")
    parser.add_argument("--keep-fastq", action="store_true", help="Keep FASTQ files after processing")
    
    args = parser.parse_args()
    
    processor = SimpleSRAProcessor(args)
    success = processor.process_all_samples()
    
    if success:
        logger.info("SRA processing completed successfully")
        return 0
    else:
        logger.error("SRA processing failed")
        return 1
        
if __name__ == "__main__":
    sys.exit(main())
