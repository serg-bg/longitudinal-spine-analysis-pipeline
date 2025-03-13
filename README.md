# Longitudinal Spine Analysis Pipeline

A comprehensive pipeline for tracking dendritic spines across longitudinal imaging data.

## Overview

This pipeline integrates spine registration, tracking, and biological classification to analyze dendritic spine dynamics across multiple timepoints.

## Features

- Robust image registration using StackReg algorithm with multiple passes
- Spine tracking based on overlap and distance metrics
- Biological classification of spines (stable, formed, eliminated, recurrent)
- Day information preservation throughout the workflow
- Batch processing support for multiple animals and genotypes
- Comprehensive visualization and statistics generation

## Installation

See the [Installation Guide](docs/installation.md) for detailed instructions.

## Usage

```bash
# Run a test on a sample segment
./run.py test

# Run batch processing on all segments
./run.py batch

# Perform a dry run to count segments without processing
./run.py dry-run
```

## Directory Structure

- `spine_tracker/`: Core Python modules
- `scripts/`: Executable scripts
- `docs/`: Documentation
- `output/`: Output directory (created on demand)

## License

MIT License

