#!/bin/bash
#SBATCH --job-name=sandfish_eval
#SBATCH --account=cs249
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/ibex/user/ahmedo/evaluation_%j.out
#SBATCH --error=/ibex/user/ahmedo/evaluation_%j.err

# Error handling function
error_exit() {
    echo "ERROR: $1" >&2
    echo "Job failed at $(date)" >&2
    exit 1
}

# Function to check if previous step succeeded
check_success() {
    if [ $? -ne 0 ]; then
        error_exit "Previous command failed"
    fi
}

echo "[$(date)] Starting assembly evaluation"

# Set up paths
ASSEMBLY_DIR="/ibex/user/ahmedo/backup_assembly/verkko"
EVAL_DIR="/ibex/user/ahmedo/evaluation"
mkdir -p ${EVAL_DIR} || error_exit "Failed to create evaluation directory"

# Input files
ASSEMBLY="${ASSEMBLY_DIR}/assembly.fasta"
HIFI="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz"
ONT="/ibex/reference/course/cs249/lizard/input/ont/lizard_ont.fastq.gz"

# Check if assembly file exists
if [ ! -f ${ASSEMBLY} ]; then
    error_exit "Assembly file not found: ${ASSEMBLY}"
fi

echo "Using assembly file: ${ASSEMBLY}"
echo "Size: $(ls -lh ${ASSEMBLY} | awk '{print $5}')"
echo "[$(date)] Running QUAST for basic metrics"
mkdir -p ${EVAL_DIR}/quast

# Load QUAST module
module purge
module load quast || error_exit "Failed to load QUAST module"

# Run QUAST
quast.py \
    -o ${EVAL_DIR}/quast \
    -t ${SLURM_CPUS_PER_TASK} \
    --large \
    ${ASSEMBLY}
check_success

echo "[$(date)] QUAST analysis completed"
echo "QUAST results:"
cat ${EVAL_DIR}/quast/report.txt