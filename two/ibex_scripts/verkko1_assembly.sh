#!/bin/bash
#SBATCH --job-name=sand_fish_assembly
#SBATCH --account=cs249
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=/ibex/user/ahmedo/assembly_%j.out
#SBATCH --error=/ibex/user/ahmedo/assembly_%j.err

# Define input files
HIFI=/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz
ONT=/ibex/reference/course/cs249/lizard/input/ont/lizard_ont.fastq.gz
HIC1=/ibex/reference/course/cs249/lizard/input/hic/lizard_hic_R1.fastq.gz
HIC2=/ibex/reference/course/cs249/lizard/input/hic/lizard_hic_R2.fastq.gz

# Create output directory
mkdir -p /ibex/user/ahmedo/lizard_assembly

# Load Verkko module
module purge
module load verkko/1.4.1

echo "[$(date)] Running Verkko"
verkko \
  -d /ibex/user/ahmedo/lizard_assembly/verkko \
  --hifi "$HIFI" \
  --nano "$ONT" \
  --hic1 "$HIC1" \
  --hic2 "$HIC2" \
  --threads 32

echo "[$(date)] Verkko assembly completed"
