#!/bin/bash
# Run batch processing for spine tracking analysis

# Determine script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Activate the virtual environment
source "${BASE_DIR}/../.venv/bin/activate"

# Default values
DATA_DIR="${BASE_DIR}/../time_tracking"
OUTPUT_DIR="${BASE_DIR}/output/batch"
GENOTYPES="B6"
ANIMALS="61-B6-027"
STACKREG_PASSES=30
QUALITY_THRESHOLD=0.85

# Process command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --genotypes)
            GENOTYPES="$2"
            shift 2
            ;;
        --animals)
            ANIMALS="$2"
            shift 2
            ;;
        --stackreg-passes)
            STACKREG_PASSES="$2"
            shift 2
            ;;
        --quality-threshold)
            QUALITY_THRESHOLD="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the batch processing
PYTHON_SCRIPT="${SCRIPT_DIR}/process_all_segments.py"

# Verify the script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Could not find the script at $PYTHON_SCRIPT"
    exit 1
fi

# Execute the script with full path
python "$PYTHON_SCRIPT" \
    --data-dir "$DATA_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --genotypes "$GENOTYPES" \
    --animals "$ANIMALS" \
    --stackreg-passes "$STACKREG_PASSES" \
    --quality-threshold "$QUALITY_THRESHOLD"

echo "Batch processing complete. Results in $OUTPUT_DIR"