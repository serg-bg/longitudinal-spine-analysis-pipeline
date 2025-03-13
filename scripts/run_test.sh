#!/bin/bash
# Run test for the integrated pipeline

# Determine script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Activate the virtual environment
source "${BASE_DIR}/../.venv/bin/activate"

# Define paths
SOURCE_DIR="${BASE_DIR}/../time_tracking/B6/61-B6-027/Segment1"
OUTPUT_DIR="${BASE_DIR}/output/test"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the test workflow
python "${SCRIPT_DIR}/test_workflow.py" \
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