# Spine Tracking Analysis Integrated Pipeline

## Overview
This pipeline integrates the functionality from multiple codebases to create a single, unified workflow for spine tracking analysis with day information preservation.

## Pipeline Components

1. **Data Loading and MIP Creation**
   - Load 3D TIFF stacks from Segmentation_Labels
   - Extract day information from filenames (metadata handling)
   - Create Maximum Intensity Projections (MIPs)

2. **Registration**
   - Apply StackReg with multiple passes for improved alignment
   - Create registered multi-channel TIFF

3. **Alignment Quality Assessment**
   - Evaluate alignment quality of each timepoint
   - Filter out poorly aligned timepoints
   - Update metadata to reflect removed timepoints

4. **Spine Tracking**
   - Extract spine objects from each timepoint
   - Track spines across timepoints
   - Classify spine status (stable, new, eliminated, etc.)
   - Preserve day information throughout tracking

5. **Visualization and Analysis**
   - Create tracking visualizations with day information
   - Generate biological classification and metrics
   - Create literature-style figures with correct day labels
   - Export CSVs with day information

6. **Reporting**
   - Generate HTML reports with alignment and tracking information
   - Include day-specific metrics

## Implementation Plan

1. Create core modules:
   - `metadata_handler.py` - For day information extraction and management
   - `spine_registration.py` - For MIP creation and registration
   - `spine_tracker.py` - For spine tracking with day information
   - `analysis.py` - For biological classification and figure generation
   - `process_all_segments.py` - Main entry point for batch processing

2. Create example workflow script showing end-to-end process
   - Include clear documentation for each step

3. Address current issues:
   - Poor alignment quality in test run
   - Ensuring proper MIP creation and registration
   - Preserving day information throughout pipeline