#!/bin/bash
#SBATCH --job-name=b_assembly
#SBATCH --account=cs249
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=/ibex/user/ahmedo/backup_assembly_%j.out
#SBATCH --error=/ibex/user/ahmedo/backup_assembly_%j.err

# Define input files
HIFI=/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz
ONT=/ibex/reference/course/cs249/lizard/input/ont/lizard_ont.fastq.gz

# Create output directory
mkdir -p /ibex/user/ahmedo/backup_assembly

# Load Verkko module
module purge
module load verkko

echo "[$(date)] Running Verkko"
verkko \
  -d /ibex/user/ahmedo/backup_assembly/verkko \
  --hifi "$HIFI" \
  --nano "$ONT" \
  --threads 32

echo "[$(date)] Verkko assembly completed"