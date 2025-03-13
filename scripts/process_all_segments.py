#!/usr/bin/env python
"""
Process All Segments Script for Spine Tracking Analysis

This script processes multiple segments in batch using the integrated pipeline.
It searches for segments in the specified directory structure and runs the 
complete workflow on each segment.
"""

import os
import glob
import argparse
import subprocess
import re
import sys

# Add the parent directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from spine_tracker.metadata_handler import find_timepoint_directories

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process multiple segments with the integrated pipeline.')
    parser.add_argument('--data-dir', '-d', default='../time_tracking', help='Path to data directory')
    parser.add_argument('--output-dir', '-o', default='./batch_output', help='Output directory')
    parser.add_argument('--genotypes', '-g', help='Comma-separated list of genotypes to process (default: all)')
    parser.add_argument('--animals', '-a', help='Comma-separated list of animals to process (default: all)')
    parser.add_argument('--segments', '-s', help='Comma-separated list of segments to process (default: all)')
    parser.add_argument('--stackreg-passes', '-p', type=int, default=30, help='Number of StackReg passes')
    parser.add_argument('--quality-threshold', '-q', type=float, default=0.85, help='Alignment quality threshold')
    parser.add_argument('--min-overlap', '-v', type=float, default=0.3, help='Minimum overlap for linking')
    parser.add_argument('--max-distance', '-x', type=float, default=10.0, help='Maximum distance for linking')
    parser.add_argument('--min-size', '-m', type=int, default=3, help='Minimum spine size in pixels')
    parser.add_argument('--skip-filter', '-f', action='store_true', help='Skip alignment quality filtering')
    parser.add_argument('--spatial-threshold', '-t', type=float, default=1.15, help='Maximum displacement in microns')
    parser.add_argument('--recurrence-window', '-r', type=int, default=3, help='Max frames a spine can disappear')
    return parser.parse_args()

def main():
    """Main function for batch processing."""
    args = parse_arguments()
    
    # Convert comma-separated lists to actual lists
    genotypes = args.genotypes.split(',') if args.genotypes else None
    animals = args.animals.split(',') if args.animals else None
    segments = args.segments.split(',') if args.segments else None
    
    # Get all timepoint directories
    timepoint_dirs = find_timepoint_directories(args.data_dir)
    
    # Filter directories based on genotype, animal, and segment
    filtered_dirs = []
    for timepoint_dir in timepoint_dirs:
        # Extract components from path
        path_parts = timepoint_dir.split(os.sep)
        
        # Try to extract genotype, animal, and segment
        genotype_match = None
        animal_match = None
        segment_match = None
        
        # Special case for SRGAP2A_Const and Fezf2_Pol2Het genotypes
        for genotype_name in ['SRGAP2A_Const', 'Fezf2_Pol2Het']:
            if genotype_name in timepoint_dir:
                genotype_match = genotype_name
                # If we found a genotype, extract the animal ID more flexibly
                for part in path_parts:
                    if genotype_name in timepoint_dir and (
                        re.match(r'SRGAP2A-\d+-\d+(-\d+)?', part) or 
                        re.match(r'FezfF2_Pol2FF-\d+', part)
                    ):
                        animal_match = part
                
        # Standard pattern matching for other genotypes
        for part in path_parts:
            # Match animal ID pattern (e.g., 61-B6-027)
            if re.match(r'\d+-[A-Za-z0-9]+-\d+', part):
                animal_match = part
            # Match segment pattern
            elif re.match(r'Segment\s*\d+', part, re.IGNORECASE) or re.match(r'apical\d+', part, re.IGNORECASE):
                segment_match = part
        
        # For genotype, check the parent directory of animal if not already set
        if not genotype_match and animal_match in path_parts:
            animal_idx = path_parts.index(animal_match)
            if animal_idx > 0:
                genotype_match = path_parts[animal_idx - 1]
        
        # Apply filters
        if genotypes and (not genotype_match or genotype_match not in genotypes):
            continue
        
        if animals and (not animal_match or not any(animal in animal_match for animal in animals)):
            continue
        
        if segments and (not segment_match or not any(segment in segment_match for segment in segments)):
            continue
        
        filtered_dirs.append(timepoint_dir)
    
    # Print summary
    print(f"Found {len(filtered_dirs)} segments to process")
    for i, directory in enumerate(filtered_dirs):
        print(f"  {i+1}. {directory}")
    
    # Process each segment
    for i, directory in enumerate(filtered_dirs):
        print(f"\nProcessing segment {i+1}/{len(filtered_dirs)}: {directory}")
        
        # Create output directory based on directory structure
        rel_path = os.path.relpath(directory, args.data_dir)
        output_dir = os.path.join(args.output_dir, rel_path)
        os.makedirs(output_dir, exist_ok=True)
        
        # Build command for test_workflow.py
        script_dir = os.path.dirname(os.path.abspath(__file__))
        test_workflow_path = os.path.join(script_dir, "test_workflow.py")
        cmd = [
            "python", test_workflow_path,
            "--source-dir", directory,
            "--output-dir", output_dir,
            "--stackreg-passes", str(args.stackreg_passes),
            "--quality-threshold", str(args.quality_threshold),
            "--min-overlap", str(args.min_overlap),
            "--max-distance", str(args.max_distance),
            "--min-size", str(args.min_size),
            "--spatial-threshold", str(args.spatial_threshold),
            "--recurrence-window", str(args.recurrence_window)
        ]
        
        if args.skip_filter:
            cmd.append("--skip-filter")
        
        # Run the workflow
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        # Print output
        print(stdout)
        if stderr:
            print("Errors:")
            print(stderr)
        
        # Check if successful
        if process.returncode == 0:
            print(f"Successfully processed: {directory}")
        else:
            print(f"Error processing: {directory}")
    
    print("\nAll segments processed!")

if __name__ == "__main__":
    main()