# Spine Tracking Analysis Pipeline

A comprehensive pipeline for tracking dendritic spines across longitudinal imaging data. This software analyzes 2D segmentation labels from microscopy data to track dendritic spines across multiple timepoints, providing quantitative analysis of spine dynamics.

## Structure
- `spine_tracker/`: Core Python modules for spine tracking and analysis
- `scripts/`: Executable scripts for running the pipeline
- `tests/`: Test cases and validation
- `docs/`: Documentation
- `output/`: Directory for output files (created on demand)

## Installation

### Using UV (Recommended)

This project uses UV (v0.5.24+), a fast Python package installer and resolver. To install:

1. Install UV if you haven't already:
   ```bash
   curl -sSf https://install.ultraviolet.rs | sh
   ```

2. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/spine-tracker.git
   cd spine-tracker
   ```

3. Run the install script (recommended):
   ```bash
   ./install.sh
   ```
   
   Or manually set up using UV commands:
   ```bash
   # Create a lockfile
   uv lock
   
   # Create a virtual environment
   uv venv
   
   # Install dependencies
   uv sync
   ```

4. To run commands with UV:
   ```bash
   # Activate the environment
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   
   # Or run directly with UV
   uv run -- ./run.py test
   ```

#### Development Dependencies

To install development dependencies for testing, linting, or documentation:

```bash
# Install test dependencies
uv pip install -e ".[test]"

# Install linting tools
uv pip install -e ".[lint]"

# Install documentation tools
uv pip install -e ".[docs]"
```

### Using pip

Alternatively, you can use pip:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/spine-tracker.git
   cd spine-tracker
   ```

2. Create a virtual environment and install dependencies:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   pip install -e .
   ```

## Usage

The pipeline can be run using the main `run.py` script:

```bash
# Run a test on a sample segment
./run.py test

# Run batch processing on all segments
./run.py batch

# Perform a dry run to count segments without processing
./run.py dry-run
```

For advanced usage and customization options, run:
```bash
./run.py --help
```

## Data Structure

The pipeline expects data to be organized in the following structure:
```
time_tracking/
├── [genotype]/               # e.g., B6, Pol2Het, etc.
│   ├── [animal-id]/          # e.g., 61-B6-027
│   │   ├── [segment-folder]/ # e.g., Segment1, apical1, etc.
│   │   │   ├── [images]      # Segmentation images with Day pattern in filename
```

## Documentation

For detailed documentation, see [docs/README.md](docs/README.md).

