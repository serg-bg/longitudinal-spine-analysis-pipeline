#!/bin/bash
# Dry run to count segments without processing

# Determine script directory - handle symlinks correctly
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

echo "Script directory: $SCRIPT_DIR"

# Activate the virtual environment
source "${BASE_DIR}/../.venv/bin/activate"

# Run the dry run for all genotypes
echo "Running dry run for all genotypes..."
PYTHON_SCRIPT="${SCRIPT_DIR}/dry_run.py"
echo "Looking for Python script at: $PYTHON_SCRIPT"

# Verify the script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Could not find the script at $PYTHON_SCRIPT"
    # Try a fallback to the scripts directory directly
    FALLBACK_SCRIPT="${BASE_DIR}/scripts/dry_run.py"
    echo "Trying fallback path: $FALLBACK_SCRIPT"
    
    if [ -f "$FALLBACK_SCRIPT" ]; then
        echo "Found script at fallback location"
        PYTHON_SCRIPT="$FALLBACK_SCRIPT"
    else
        echo "Error: Could not find script at fallback location either"
        exit 1
    fi
fi

python "$PYTHON_SCRIPT" --data-dir "${BASE_DIR}/../time_tracking" --min-timepoints 3

# Run the dry run for specific genotypes
echo -e "\nRunning dry run for B6 genotype..."
python "$PYTHON_SCRIPT" --data-dir "${BASE_DIR}/../time_tracking" --genotypes B6 --min-timepoints 3

echo -e "\nRunning dry run for Pol2Het genotype..."
python "$PYTHON_SCRIPT" --data-dir "${BASE_DIR}/../time_tracking" --genotypes Pol2Het --min-timepoints 3

echo -e "\nRunning dry run for Pol3Het genotype..."
python "$PYTHON_SCRIPT" --data-dir "${BASE_DIR}/../time_tracking" --genotypes Pol3Het --min-timepoints 3

echo -e "\nRunning dry run for SRGAP2A_Const genotype..."
python "$PYTHON_SCRIPT" --data-dir "${BASE_DIR}/../time_tracking" --genotypes SRGAP2A_Const --min-timepoints 3

echo -e "\nRunning dry run for Fezf2_Pol2Het genotype..."
python "$PYTHON_SCRIPT" --data-dir "${BASE_DIR}/../time_tracking" --genotypes Fezf2_Pol2Het --min-timepoints 3

echo "Dry run completed."