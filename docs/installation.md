# Installation Guide

## Prerequisites

- Python 3.11 or later
- UV package manager (recommended) or pip

## Using UV (Recommended)

1. Install UV if you haven't already:
   ```bash
   curl -sSf https://install.ultraviolet.rs | sh
   ```

2. Clone the repository:
   ```bash
   git clone https://github.com/sergiobgarcia/longitudinal-spine-analysis-pipeline.git
   cd longitudinal-spine-analysis-pipeline
   ```

3. Run the install script:
   ```bash
   ./install.sh
   ```

## Using pip

1. Clone the repository:
   ```bash
   git clone https://github.com/sergiobgarcia/longitudinal-spine-analysis-pipeline.git
   cd longitudinal-spine-analysis-pipeline
   ```

2. Create a virtual environment and install dependencies:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scriptsctivate
   pip install -e .
   ```

## Development Dependencies

To install development dependencies for testing, linting, or documentation:

```bash
# Install test dependencies
pip install -e ".[test]"

# Install linting tools
pip install -e ".[lint]"

# Install documentation tools
pip install -e ".[docs]"
```

