# Spine Tracking Analysis Integrated Pipeline

This pipeline provides a unified workflow for tracking dendritic spines across timepoints in longitudinal imaging data. It combines functionality from multiple existing codebases and adds robust day information preservation throughout the analysis.

## Features

- **Metadata Extraction**: Extract day information from file naming schemes
- **MIP Creation**: Generate Maximum Intensity Projections (MIPs) from 3D TIFF stacks
- **Registration**: Apply StackReg with multiple passes to improve alignment
- **Alignment Quality Assessment**: Evaluate and filter timepoints based on alignment quality
- **Spine Tracking**: Track spines across timepoints with day information preservation
- **Visualization**: Generate visualizations with actual experimental day labels

## Pipeline Components

The pipeline consists of the following key modules:

1. **metadata_handler.py**: Handles extraction and management of metadata
2. **spine_registration.py**: Creates MIPs and handles registration
3. **spine_tracker.py**: Tracks spines and preserves day information
4. **test_workflow.py**: Demonstrates basic pipeline usage

## Quick Start

To run a test of the pipeline:

```bash
# Clone the repository
git clone https://github.com/your-repo/spine-tracking-analysis.git
cd spine-tracking-analysis/integrated_pipeline

# Install dependencies
pip install -r requirements.txt

# Run the test workflow
./run_test.sh
```

## Usage

The general workflow is:

1. Extract metadata from segmentation files
2. Create MIPs and register them using StackReg
3. Assess alignment quality and filter poor timepoints
4. Track spines across timepoints
5. Generate visualizations and analysis figures

### Example

```python
# Import pipeline modules
from metadata_handler import MetadataHandler
import spine_registration

# Initialize metadata handler
metadata_handler = MetadataHandler()
metadata_handler.extract_from_segmentation_files(source_dir)

# Process the segment
segmentation_path, good_frames = spine_registration.process_segment(
    source_dir=source_dir,
    output_dir=output_dir,
    num_stackreg_passes=30,
    quality_threshold=0.85
)

# Update metadata if frames were filtered
if len(good_frames) < len(metadata_handler.metadata['timepoints']):
    metadata_handler.update_removed_frames(good_frames)

# Save metadata
metadata_path = metadata_handler.save_metadata(output_dir)
```

## Requirements

- Python 3.7+
- NumPy
- pandas
- tifffile
- scikit-image
- matplotlib
- pystackreg (for StackReg functionality)

## License

MIT License

## Acknowledgments

This pipeline builds upon work from multiple sources and teams.