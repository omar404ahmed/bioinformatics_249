#!/bin/bash
#SBATCH --job-name=inspector
#SBATCH --account=cs249
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/ibex/user/ahmedo/inspector_%j.out
#SBATCH --error=/ibex/user/ahmedo/inspector_%j.err

echo "[$(date)] Starting Inspector evaluation"

# Set paths
ASSEMBLY_DIR="/ibex/user/ahmedo/backup_assembly/verkko"
EVAL_DIR="/ibex/user/ahmedo/evaluation/inspector"
mkdir -p ${EVAL_DIR}

# Input files
ASSEMBLY="${ASSEMBLY_DIR}/assembly.fasta"
HIFI="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz"

# Change to working directory
cd ${EVAL_DIR}

# Load Inspector module
module purge
module load inspector

echo "[$(date)] Running Inspector to detect mis-assemblies"

# Run Inspector
inspector.py \
    -c ${ASSEMBLY} \
    -r ${HIFI} \
    -o ${EVAL_DIR}/output \
    -t ${SLURM_CPUS_PER_TASK}

# Check if Inspector completed successfully
if [ $? -eq 0 ]; then
    echo "[$(date)] Inspector analysis completed successfully"
    
    # Display summary if available
    if [ -f "${EVAL_DIR}/output/summary.txt" ]; then
        echo "Inspector results summary:"
        cat ${EVAL_DIR}/output/summary.txt
    fi
else
    echo "[$(date)] Inspector analysis failed"
fi

echo "[$(date)] Inspector evaluation complete"