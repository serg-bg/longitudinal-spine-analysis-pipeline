#!/usr/bin/env python
"""
Test Workflow for Spine Tracking Analysis Integrated Pipeline

This script demonstrates a complete workflow using the integrated pipeline
to process a single segment from start to finish, including:
1. Metadata extraction
2. MIP creation and registration
3. Alignment quality assessment
4. Spine tracking
5. Biological classification
6. Figure generation
"""

import os
import sys
import argparse

# Add the parent directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from spine_tracker.metadata_handler import MetadataHandler
from spine_tracker import spine_registration
from spine_tracker import spine_tracking as spine_tracker
from spine_tracker import analysis

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Test workflow for spine tracking analysis pipeline.')
    parser.add_argument('--source-dir', required=True, help='Source directory with Segmentation_Labels')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--stackreg-passes', type=int, default=30, help='Number of StackReg passes')
    parser.add_argument('--quality-threshold', type=float, default=0.85, help='Alignment quality threshold')
    parser.add_argument('--min-overlap', type=float, default=0.3, help='Minimum overlap for linking spines')
    parser.add_argument('--max-distance', type=float, default=10.0, help='Maximum distance for linking spines')
    parser.add_argument('--min-size', type=int, default=3, help='Minimum spine size in pixels')
    parser.add_argument('--skip-filter', action='store_true', help='Skip alignment quality filtering')
    parser.add_argument('--spatial-threshold', type=float, default=1.15, help='Maximum spine displacement in microns')
    parser.add_argument('--recurrence-window', type=int, default=3, help='Max frames a spine can disappear')
    return parser.parse_args()

def main():
    """Main function for test workflow."""
    args = parse_arguments()
    
    print("=== Spine Tracking Analysis Integrated Pipeline ===")
    print(f"Source directory: {args.source_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"StackReg passes: {args.stackreg_passes}")
    print(f"Quality threshold: {args.quality_threshold}")
    print(f"Skip filter: {args.skip_filter}")
    print(f"Tracking parameters:")
    print(f"  - Min spine size: {args.min_size} pixels")
    print(f"  - Max distance: {args.max_distance} pixels")
    print(f"  - Min overlap: {args.min_overlap}")
    print(f"Analysis parameters:")
    print(f"  - Spatial threshold: {args.spatial_threshold} microns")
    print(f"  - Recurrence window: {args.recurrence_window} frames")
    print("===================================================")
    
    # Step 1: Initialize metadata handler and extract metadata
    print("\nStep 1: Extracting metadata...")
    metadata_handler = MetadataHandler()
    metadata_handler.extract_from_path(args.source_dir)
    metadata_handler.extract_from_segmentation_files(args.source_dir)
    
    print(f"  Animal ID: {metadata_handler.metadata['animal_id']}")
    print(f"  Genotype: {metadata_handler.metadata['genotype']}")
    print(f"  Segment: {metadata_handler.metadata['segment']}")
    
    if metadata_handler.metadata['timepoints']:
        print(f"  Timepoints: {len(metadata_handler.metadata['timepoints'])}")
        for frame, info in metadata_handler.metadata['timepoints'].items():
            print(f"    Frame {frame} = Day {info['day']}")
    else:
        print("  No day information found.")
    
    # Step 2: Process the segment to create registered segmentation
    print("\nStep 2: Creating registered segmentation...")
    segmentation_path, good_frames = spine_registration.process_segment(
        source_dir=args.source_dir,
        output_dir=args.output_dir,
        num_stackreg_passes=args.stackreg_passes,
        quality_threshold=args.quality_threshold,
        skip_filter=args.skip_filter
    )
    
    # Step 3: Update metadata if frames were filtered
    if len(good_frames) < len(metadata_handler.metadata['timepoints']):
        print("\nStep 3: Updating metadata for filtered frames...")
        metadata_handler.update_removed_frames(good_frames)
        
        print(f"  Removed frames: {metadata_handler.metadata['removed_frames']}")
        print(f"  Removed days: {metadata_handler.metadata['removed_days']}")
        
        # Updated timepoints
        print(f"  Updated timepoints: {len(metadata_handler.metadata['timepoints'])}")
        for frame, info in metadata_handler.metadata['timepoints'].items():
            print(f"    Frame {frame} = Day {info['day']}")
    
    # Step 4: Save metadata
    print("\nStep 4: Saving metadata...")
    metadata_path = metadata_handler.save_metadata(args.output_dir)
    
    # Step 5: Track spines
    print("\nStep 5: Tracking spines...")
    spine_features, track_labels = spine_tracker.track_spines_with_metadata(
        segmentation_path=segmentation_path,
        metadata=metadata_handler.metadata,
        output_dir=args.output_dir,
        min_size=args.min_size,
        max_distance=args.max_distance,
        min_overlap=args.min_overlap
    )
    
    # Step 6: Apply biological criteria and create visualizations
    print("\nStep 6: Applying biological criteria and creating visualizations...")
    validation_dir = os.path.join(args.output_dir, 'validation')
    analysis_results = analysis.process_tracking_results(
        input_dir=args.output_dir,
        output_dir=validation_dir,
        spatial_threshold=args.spatial_threshold,
        recurrence_window=args.recurrence_window
    )
    
    print("\nWorkflow completed successfully!")
    print(f"Registered segmentation: {segmentation_path}")
    print(f"Spine tracks: {os.path.join(args.output_dir, os.path.basename(segmentation_path).replace('.tif', '_spine_tracks.csv'))}")
    print(f"Metadata file: {metadata_path}")
    print(f"Validation report: {os.path.join(validation_dir, 'biological_criteria_report.html')}")
    
if __name__ == '__main__':
    main()