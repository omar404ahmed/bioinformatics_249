#!/bin/bash
#SBATCH --job-name=merqury_eval
#SBATCH --account=cs249
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/ibex/user/ahmedo/merqury_%j.out
#SBATCH --error=/ibex/user/ahmedo/merqury_%j.err

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

echo "[$(date)] Starting Merqury evaluation"

# Set up paths
ASSEMBLY_DIR="/ibex/user/ahmedo/backup_assembly/verkko"
EVAL_DIR="/ibex/user/ahmedo/merqury"
mkdir -p ${EVAL_DIR} || error_exit "Failed to create evaluation directory"

# Input files
ASSEMBLY="${ASSEMBLY_DIR}/assembly.fasta"
HIFI="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz"

# Check if files exist
if [ ! -f ${ASSEMBLY} ]; then
    error_exit "Assembly file not found: ${ASSEMBLY}"
fi
if [ ! -f ${HIFI} ]; then
    error_exit "HiFi reads file not found: ${HIFI}"
fi

echo "Using assembly file: ${ASSEMBLY}"
echo "Using HiFi reads: ${HIFI}"

cd ${EVAL_DIR} || error_exit "Failed to change to evaluation directory"

# Check if meryl/merqury is available as a module, otherwise use conda
if module load meryl merqury &>/dev/null; then
    echo "Using meryl and merqury from modules"
    MERYL_CMD="meryl"
    MERQURY_CMD="merqury.sh"
else
    echo "Meryl/Merqury modules not found, setting up conda environment"
    
    # Check if conda is available
    if ! command -v conda &>/dev/null; then
        error_exit "Neither modules nor conda is available for meryl/merqury"
    fi
    
    # Set up conda environment
    CONDA_ENV="${EVAL_DIR}/conda_env"
    conda create -p ${CONDA_ENV} -c bioconda -c conda-forge merqury python=3.8 -y || error_exit "Failed to create conda environment"
    
    # Activate conda environment
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate ${CONDA_ENV} || error_exit "Failed to activate conda environment"
    
    MERYL_CMD="meryl"
    MERQURY_CMD="merqury.sh"
    # Verify commands are available
    if ! command -v ${MERYL_CMD} &>/dev/null; then
        error_exit "meryl command not found after setting up conda"
    fi
    if ! command -v ${MERQURY_CMD} &>/dev/null; then
        error_exit "merqury.sh command not found after setting up conda"
    fi
fi

# Create output directory for k-mer database
MERYL_DB="${EVAL_DIR}/meryl_db"
mkdir -p ${MERYL_DB} || error_exit "Failed to create meryl database directory"

# Run meryl to build k-mer database
echo "[$(date)] Building k-mer database with meryl"
K=21  # Default k-mer size, can be adjusted
${MERYL_CMD} count k=${K} ${HIFI} output ${MERYL_DB} || error_exit "Failed to run meryl count"
check_success

# Print k-mer statistics
echo "[$(date)] K-mer statistics:"
${MERYL_CMD} statistics ${MERYL_DB} > ${EVAL_DIR}/kmer_statistics.txt
check_success
cat ${EVAL_DIR}/kmer_statistics.txt

# Run merqury
echo "[$(date)] Running merqury for QV evaluation"
OUT_PREFIX="${EVAL_DIR}/sandfish_merqury"
${MERQURY_CMD} ${MERYL_DB} ${ASSEMBLY} ${OUT_PREFIX} || error_exit "Failed to run merqury"
check_success

echo "[$(date)] Merqury analysis completed"

# Extract and display QV score
if [ -f "${OUT_PREFIX}.qv" ]; then
    echo "QV score results:"
    cat ${OUT_PREFIX}.qv
    
    # Extract key metrics for the report
    QV_SCORE=$(grep "QV" ${OUT_PREFIX}.qv | awk '{print $2}')
    ERROR_RATE=$(grep "ERROR" ${OUT_PREFIX}.qv | awk '{print $2}')
else
    echo "WARNING: QV score file not found at ${OUT_PREFIX}.qv"
    QV_SCORE="N/A"
    ERROR_RATE="N/A"
fi

# Create summary report
echo "[$(date)] Creating Merqury summary report"

cat > ${EVAL_DIR}/merqury_summary.txt << EOF
==============================================================
Scincus mitranus (Sandfish) Genome Assembly - Merqury Evaluation
==============================================================
Date: $(date)

K-MER ANALYSIS AND QV SCORE
---------------------------
Assembly file: ${ASSEMBLY}
HiFi reads: ${HIFI}
K-mer size: ${K}

QV Score: ${QV_SCORE}
(QV is a log-scaled measure of the assembly's base-level accuracy)

Estimated error rate: ${ERROR_RATE}
(Lower error rate indicates higher accuracy)

EOF

echo "[$(date)] Merqury evaluation completed. Summary report available at: ${EVAL_DIR}/merqury_summary.txt"