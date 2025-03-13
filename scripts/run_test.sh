#!/bin/bash
# Run test for the integrated pipeline

# Determine script directory - handle symlinks correctly
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

echo "Script directory: $SCRIPT_DIR"

# Activate the virtual environment
source "${BASE_DIR}/../.venv/bin/activate"

# Define paths
SOURCE_DIR="${BASE_DIR}/../time_tracking/B6/61-B6-027/Segment1"
OUTPUT_DIR="${BASE_DIR}/output/test"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the test workflow
PYTHON_SCRIPT="${SCRIPT_DIR}/test_workflow.py"
echo "Looking for Python script at: $PYTHON_SCRIPT"

# Verify the script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Could not find the script at $PYTHON_SCRIPT"
    # Try a fallback to the scripts directory directly
    FALLBACK_SCRIPT="${BASE_DIR}/scripts/test_workflow.py"
    echo "Trying fallback path: $FALLBACK_SCRIPT"
    
    if [ -f "$FALLBACK_SCRIPT" ]; then
        echo "Found script at fallback location"
        PYTHON_SCRIPT="$FALLBACK_SCRIPT"
    else
        echo "Error: Could not find script at fallback location either"
        exit 1
    fi
fi

python "$PYTHON_SCRIPT" \
    --source-dir "$SOURCE_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --stackreg-passes 30 \
    --quality-threshold 0.85 \
    --min-overlap 0.3 \
    --max-distance 10.0 \
    --min-size 3 \
    --spatial-threshold 1.15 \
    --recurrence-window 3

echo "Test complete. Results in $OUTPUT_DIR"
echo "HTML report available at: $OUTPUT_DIR/validation/biological_criteria_report.html"