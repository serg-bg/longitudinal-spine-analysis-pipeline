#!/bin/bash
# Dry run to count segments without processing

# Determine script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Activate the virtual environment
source "${BASE_DIR}/../.venv/bin/activate"

# Run the dry run for all genotypes
echo "Running dry run for all genotypes..."
python "${SCRIPT_DIR}/dry_run.py" --data-dir "${BASE_DIR}/../time_tracking" --min-timepoints 3

# Run the dry run for specific genotypes
echo -e "\nRunning dry run for B6 genotype..."
python "${SCRIPT_DIR}/dry_run.py" --data-dir "${BASE_DIR}/../time_tracking" --genotypes B6 --min-timepoints 3

echo -e "\nRunning dry run for Pol2Het genotype..."
python "${SCRIPT_DIR}/dry_run.py" --data-dir "${BASE_DIR}/../time_tracking" --genotypes Pol2Het --min-timepoints 3

echo -e "\nRunning dry run for Pol3Het genotype..."
python "${SCRIPT_DIR}/dry_run.py" --data-dir "${BASE_DIR}/../time_tracking" --genotypes Pol3Het --min-timepoints 3

echo -e "\nRunning dry run for SRGAP2A_Const genotype..."
python "${SCRIPT_DIR}/dry_run.py" --data-dir "${BASE_DIR}/../time_tracking" --genotypes SRGAP2A_Const --min-timepoints 3

echo -e "\nRunning dry run for Fezf2_Pol2Het genotype..."
python "${SCRIPT_DIR}/dry_run.py" --data-dir "${BASE_DIR}/../time_tracking" --genotypes Fezf2_Pol2Het --min-timepoints 3

echo "Dry run completed."