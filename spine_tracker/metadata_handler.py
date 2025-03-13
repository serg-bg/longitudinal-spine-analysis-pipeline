#!/usr/bin/env python
"""
Metadata Handler for Spine Tracking Analysis

This module handles the extraction and management of metadata from filenames
and directory structures in the spine tracking pipeline, ensuring that day information
and experimental metadata is preserved throughout the analysis process.
"""

import os
import re
import pandas as pd
import json

# Regular expression patterns for extracting metadata
ANIMAL_ID_PATTERN = r'(\d+-[A-Za-z0-9]+-\d+)'  # e.g., 61-B6-027
DAY_PATTERN = r'Day(\d+)'  # e.g., Day2
APICAL_PATTERN = r'apical(\d+)'  # e.g., apical1
SEGMENT_PATTERN = r'Segment\s*(\d+)'  # e.g., Segment1 or Segment 1
GENOTYPE_PATTERN = r'([A-Za-z0-9_]+)(?=/\d+-)'  # Extracts genotype from directory structure

class MetadataHandler:
    """
    Handles extraction and management of metadata from spine tracking files.
    """
    
    def __init__(self):
        """Initialize metadata handler."""
        self.metadata = {
            'animal_id': None,
            'genotype': None,
            'segment': None,
            'timepoints': {},  # Maps frame_index to day number
            'day_to_frame': {},  # Maps day number to frame_index
            'frame_to_day': {},  # Maps frame_index to day number
            'removed_frames': [],  # List of frames removed during processing
            'removed_days': [],  # List of days removed during processing
            'microns_per_pixel': 0.102,  # Default scale
        }
    
    def extract_from_path(self, filepath):
        """
        Extract metadata from a file path.
        
        Args:
            filepath: Path to a file or directory
            
        Returns:
            Dictionary of extracted metadata
        """
        # Convert to absolute path if not already
        filepath = os.path.abspath(filepath)
        
        # Extract animal ID
        animal_id_match = re.search(ANIMAL_ID_PATTERN, filepath)
        if animal_id_match:
            self.metadata['animal_id'] = animal_id_match.group(1)
            
        # Extract genotype from directory structure
        genotype_match = re.search(GENOTYPE_PATTERN, filepath)
        if genotype_match:
            self.metadata['genotype'] = genotype_match.group(1)
            
        # Extract segment number
        segment_match = re.search(SEGMENT_PATTERN, filepath)
        if segment_match:
            self.metadata['segment'] = f"Segment{segment_match.group(1)}"
            
        # Extract apical dendrite number (alternative to segment)
        apical_match = re.search(APICAL_PATTERN, filepath)
        if apical_match and not segment_match:
            self.metadata['segment'] = f"Apical{apical_match.group(1)}"
            
        return self.metadata
    
    def extract_from_segmentation_files(self, directory):
        """
        Extract day information from segmentation label files in a directory.
        
        Args:
            directory: Path to directory containing segmentation files
            
        Returns:
            Updated metadata dictionary
        """
        # Find Segmentation_Labels directory
        seg_labels_dir = None
        if os.path.basename(directory) == 'Segmentation_Labels':
            seg_labels_dir = directory
        else:
            potential_dir = os.path.join(directory, 'Segmentation_Labels')
            if os.path.exists(potential_dir):
                seg_labels_dir = potential_dir
        
        if not seg_labels_dir:
            print("Warning: Could not find Segmentation_Labels directory")
            return self.metadata
        
        # Get all TIFF files in the directory
        tiff_files = [f for f in os.listdir(seg_labels_dir) if f.endswith('.tif') or f.endswith('.tiff')]
        
        # Extract day information from each file
        day_info = []
        for tiff_file in tiff_files:
            day_match = re.search(DAY_PATTERN, tiff_file)
            if day_match:
                day_num = int(day_match.group(1))
                day_info.append((tiff_file, day_num))
                
                # Also extract animal ID if not already done
                if not self.metadata['animal_id']:
                    animal_id_match = re.search(ANIMAL_ID_PATTERN, tiff_file)
                    if animal_id_match:
                        self.metadata['animal_id'] = animal_id_match.group(1)
        
        # Sort by day number
        day_info.sort(key=lambda x: x[1])
        
        # Create mapping between frame index and day number
        for i, (file, day) in enumerate(day_info):
            self.metadata['timepoints'][i] = {
                'day': day,
                'filename': file
            }
            self.metadata['day_to_frame'][day] = i
            self.metadata['frame_to_day'][i] = day
        
        return self.metadata
    
    def update_removed_frames(self, good_frames, original_frames=None):
        """
        Update metadata to reflect frames that were removed during alignment filtering.
        
        Args:
            good_frames: List of frame indices that were kept
            original_frames: Optional list of all original frames (defaults to range(max_frame+1))
            
        Returns:
            Updated metadata dictionary
        """
        if original_frames is None:
            if self.metadata['timepoints']:
                original_frames = list(range(max(self.metadata['frame_to_day'].keys()) + 1))
            else:
                # Cannot determine original frames
                return self.metadata
            
        # Find removed frames
        removed_frames = [f for f in original_frames if f not in good_frames]
        self.metadata['removed_frames'] = removed_frames
        
        # Update day-to-frame and frame-to-day mappings
        new_frame_to_day = {}
        new_day_to_frame = {}
        new_timepoints = {}
        
        # Create mapping from old frame indices to new frame indices
        old_to_new_frame = {old_frame: new_frame for new_frame, old_frame in enumerate(good_frames)}
        
        # Update day mappings
        for old_frame, new_frame in old_to_new_frame.items():
            if old_frame in self.metadata['frame_to_day']:
                day = self.metadata['frame_to_day'][old_frame]
                new_frame_to_day[new_frame] = day
                new_day_to_frame[day] = new_frame
                
                if old_frame in self.metadata['timepoints']:
                    new_timepoints[new_frame] = self.metadata['timepoints'][old_frame]
        
        # Update metadata with new mappings
        self.metadata['frame_to_day'] = new_frame_to_day
        self.metadata['day_to_frame'] = new_day_to_frame
        self.metadata['timepoints'] = new_timepoints
        
        # Update removed days
        all_days = set(self.metadata['day_to_frame'].keys())
        remaining_days = set(new_day_to_frame.keys())
        removed_days = all_days - remaining_days
        self.metadata['removed_days'] = sorted(list(removed_days))
        
        return self.metadata
    
    def save_metadata(self, output_dir):
        """
        Save metadata to JSON file.
        
        Args:
            output_dir: Directory to save the metadata file
            
        Returns:
            Path to the saved metadata file
        """
        os.makedirs(output_dir, exist_ok=True)
        metadata_path = os.path.join(output_dir, 'tracking_metadata.json')
        
        with open(metadata_path, 'w') as f:
            json.dump(self.metadata, f, indent=2)
            
        print(f"Saved metadata to {metadata_path}")
        return metadata_path
    
    def load_metadata(self, metadata_path):
        """
        Load metadata from JSON file.
        
        Args:
            metadata_path: Path to the metadata JSON file
            
        Returns:
            Loaded metadata dictionary
        """
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                self.metadata = json.load(f)
                
            # Convert string keys to integers where needed
            self.metadata['timepoints'] = {int(k): v for k, v in self.metadata['timepoints'].items()}
            self.metadata['day_to_frame'] = {int(k): int(v) for k, v in self.metadata['day_to_frame'].items()}
            self.metadata['frame_to_day'] = {int(k): int(v) for k, v in self.metadata['frame_to_day'].items()}
            
            print(f"Loaded metadata from {metadata_path}")
        else:
            print(f"Metadata file not found: {metadata_path}")
            
        return self.metadata
    
    def get_day_labels(self):
        """
        Get day labels for plots.
        
        Returns:
            List of day labels for each frame
        """
        if not self.metadata['frame_to_day']:
            return None
            
        max_frame = max(self.metadata['frame_to_day'].keys())
        day_labels = ["" for _ in range(max_frame + 1)]
        
        for frame, day in self.metadata['frame_to_day'].items():
            day_labels[frame] = f"Day {day}"
            
        return day_labels
    
def find_timepoint_directories(base_dir):
    """
    Find all timepoint directories in the time_tracking hierarchy.
    
    Args:
        base_dir: Base directory to start the search
        
    Returns:
        List of directories containing timepoint data
    """
    timepoint_dirs = []
    
    for root, dirs, files in os.walk(base_dir):
        if 'Segmentation_Labels' in dirs:
            # This directory contains timepoint data
            timepoint_dirs.append(root)
    
    return timepoint_dirs