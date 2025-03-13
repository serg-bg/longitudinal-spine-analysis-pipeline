#!/bin/bash
# Run batch processing for spine tracking analysis

# Determine script directory - handle symlinks correctly
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

echo "Script directory: $SCRIPT_DIR"

# Activate the virtual environment
source "${BASE_DIR}/../.venv/bin/activate"

# Default values
DATA_DIR="${BASE_DIR}/../time_tracking"
OUTPUT_DIR="${BASE_DIR}/output/batch"
GENOTYPES="B6"
ANIMALS=""  # Default to all animals
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
echo "Looking for Python script at: $PYTHON_SCRIPT"

# Verify the script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Could not find the script at $PYTHON_SCRIPT"
    # Try a fallback to the scripts directory directly
    FALLBACK_SCRIPT="${BASE_DIR}/scripts/process_all_segments.py"
    echo "Trying fallback path: $FALLBACK_SCRIPT"
    
    if [ -f "$FALLBACK_SCRIPT" ]; then
        echo "Found script at fallback location"
        PYTHON_SCRIPT="$FALLBACK_SCRIPT"
    else
        echo "Error: Could not find script at fallback location either"
        exit 1
    fi
fi

# Execute the script with full path
# Build command based on provided parameters
CMD=("python" "$PYTHON_SCRIPT" 
     "--data-dir" "$DATA_DIR" 
     "--output-dir" "$OUTPUT_DIR" 
     "--genotypes" "$GENOTYPES" 
     "--stackreg-passes" "$STACKREG_PASSES" 
     "--quality-threshold" "$QUALITY_THRESHOLD")

# Only add animals parameter if specified
if [ ! -z "$ANIMALS" ]; then
    CMD+=("--animals" "$ANIMALS")
    echo "Filtering by animals: $ANIMALS"
else
    echo "Processing all animals"
fi

echo "Running command: ${CMD[@]}"
"${CMD[@]}"

echo "Batch processing complete. Results in $OUTPUT_DIR"