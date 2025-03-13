#!/usr/bin/env python
"""
Dry run script for counting the segments that would be processed
without actually running the registration and analysis.
"""

import os
import re
import argparse
import sys

# Add the parent directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from spine_tracker.metadata_handler import find_timepoint_directories

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Count segments for processing.')
    parser.add_argument('--data-dir', '-d', default='../time_tracking', help='Path to data directory')
    parser.add_argument('--genotypes', '-g', help='Comma-separated list of genotypes to process (default: all)')
    parser.add_argument('--animals', '-a', help='Comma-separated list of animals to process (default: all)')
    parser.add_argument('--segments', '-s', help='Comma-separated list of segments to process (default: all)')
    parser.add_argument('--min-timepoints', '-t', type=int, default=3, help='Minimum timepoints required (default: 3)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    return parser.parse_args()

def main():
    """Main function for dry run."""
    args = parse_arguments()
    
    # Convert comma-separated lists to actual lists
    genotypes = args.genotypes.split(',') if args.genotypes else None
    animals = args.animals.split(',') if args.animals else None
    segments = args.segments.split(',') if args.segments else None
    
    # Convert relative path to absolute path if needed
    data_dir = os.path.abspath(args.data_dir)
    
    print(f"Dry run on data directory: {data_dir}")
    print(f"Genotypes filter: {genotypes}")
    print(f"Animals filter: {animals}")
    print(f"Segments filter: {segments}")
    print(f"Minimum timepoints: {args.min_timepoints}")
    print("\nSearching for segments...")
    
    # Find all timepoint directories
    all_timepoint_dirs = find_timepoint_directories(data_dir)
    
    if args.verbose:
        print(f"Found {len(all_timepoint_dirs)} total timepoint directories")
    
    # Filter directories based on genotype, animal, and segment
    filtered_dirs = []
    for directory in all_timepoint_dirs:
        # Extract components from path
        path_parts = directory.split(os.sep)
        
        # Try to extract genotype, animal, and segment
        genotype_match = None
        animal_match = None
        segment_match = None
        
        # Special case for SRGAP2A_Const and Fezf2_Pol2Het genotypes
        if 'SRGAP2A_Const' in directory:
            genotype_match = 'SRGAP2A_Const'
            for part in path_parts:
                if re.match(r'SRGAP2A-\d+-\d+(-\d+)?', part):
                    animal_match = part
        elif 'Fezf2_Pol2Het' in directory:
            genotype_match = 'Fezf2_Pol2Het'
            for part in path_parts:
                if re.match(r'FezfF2_Pol2FF-\d+', part):
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
        
        filtered_dirs.append(directory)
    
    # Use filtered directories for the rest of the analysis
    timepoint_dirs = filtered_dirs
    
    # Count total directories and those with sufficient timepoints
    total_dirs = len(timepoint_dirs)
    valid_dirs = 0
    genotype_counts = {}
    animal_counts = {}
    segment_counts = {}
    timepoint_distribution = {}
    
    for directory in timepoint_dirs:
        # Count timepoints (number of .tif files)
        tif_files = [f for f in os.listdir(directory) if f.endswith('.tif') and 'Day' in f]
        timepoint_count = len(tif_files)
        
        # Parse directory path to extract components
        parts = directory.split(os.sep)
        genotype = 'unknown'
        animal = 'unknown'
        segment = 'unknown'
        
        # Extract genotype, animal, and segment from path
        for i, part in enumerate(parts):
            if part in ['B6', 'Pol2Het', 'Pol3Het', 'SRGAP2A_Const', 'Fezf2_Pol2Het']:
                genotype = part
                if i + 1 < len(parts):
                    animal = parts[i + 1]
                if i + 2 < len(parts):
                    segment = parts[i + 2]
        
        # Update counts
        genotype_counts[genotype] = genotype_counts.get(genotype, 0) + 1
        animal_counts[animal] = animal_counts.get(animal, 0) + 1
        segment_counts[segment] = segment_counts.get(segment, 0) + 1
        
        # Track timepoint distribution
        timepoint_distribution[timepoint_count] = timepoint_distribution.get(timepoint_count, 0) + 1
        
        # Check if it has sufficient timepoints
        if timepoint_count >= args.min_timepoints:
            valid_dirs += 1
            print(f"  - {directory} ({timepoint_count} timepoints)")
    
    # Print summary
    print("\nSummary:")
    print(f"Total segments found: {total_dirs}")
    print(f"Segments with >= {args.min_timepoints} timepoints: {valid_dirs}")
    
    print("\nGenotype distribution:")
    for genotype, count in sorted(genotype_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  - {genotype}: {count} segments")
    
    print("\nAnimal distribution:")
    for animal, count in sorted(animal_counts.items(), key=lambda x: x[1], reverse=True)[:10]:  # Top 10
        print(f"  - {animal}: {count} segments")
    if len(animal_counts) > 10:
        print(f"  - ... and {len(animal_counts) - 10} more animals")
    
    print("\nSegment name distribution:")
    for segment, count in sorted(segment_counts.items(), key=lambda x: x[1], reverse=True)[:10]:  # Top 10
        print(f"  - {segment}: {count} instances")
    if len(segment_counts) > 10:
        print(f"  - ... and {len(segment_counts) - 10} more segment names")
    
    print("\nTimepoint distribution:")
    for timepoints, count in sorted(timepoint_distribution.items()):
        print(f"  - {timepoints} timepoints: {count} segments")

if __name__ == "__main__":
    main()