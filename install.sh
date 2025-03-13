#!/bin/bash
# Installation script for Spine Tracking Analysis Pipeline

# Ensure we're in the right directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "${SCRIPT_DIR}"

# Check if UV is installed
if ! command -v uv &> /dev/null; then
    echo "UV is not installed. Would you like to install it? (y/n)"
    read -r response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        echo "Installing UV..."
        curl -sSf https://install.ultraviolet.rs | sh
    else
        echo "UV is required for installation. Please install it manually and try again."
        echo "You can install UV with: curl -sSf https://install.ultraviolet.rs | sh"
        exit 1
    fi
fi

# Check UV version
UV_VERSION=$(uv --version | cut -d' ' -f2)
echo "Using UV version: $UV_VERSION"

# Check Python version
PYTHON_VERSION=$(python --version 2>&1 | cut -d' ' -f2)
echo "Using Python version: $PYTHON_VERSION"

# Create a lockfile
echo "Creating lockfile..."
uv lock

# Create virtual environment
echo "Creating virtual environment..."
uv venv

# Sync dependencies to the environment
echo "Syncing dependencies..."
uv sync

echo ""
echo "Installation complete! You can activate the virtual environment with:"
if [[ "$OSTYPE" == "darwin"* ]] || [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "  source .venv/bin/activate"
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
    echo "  .venv\\Scripts\\activate"
fi
echo ""
echo "To run the pipeline:"
echo "  uv run -- ./run.py --help"
echo "or activate the environment and run directly:"
echo "  ./run.py --help"