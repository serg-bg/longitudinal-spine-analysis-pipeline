#!/usr/bin/env python
"""
Spine Tracking Module for Spine Tracking Analysis

This module handles tracking of dendritic spines across timepoints with
day information preservation.
"""

import os
import numpy as np
import pandas as pd
import tifffile
from skimage import measure
from skimage.morphology import binary_dilation
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

def extract_spine_objects(segmentation, min_size=3):
    """
    Extract spine objects from a segmentation image.
    
    Args:
        segmentation: 2D segmentation image (0=background, 1=spines, 2=dendrite)
        min_size: Minimum spine size in pixels
        
    Returns:
        labeled_spines: 2D labeled image of spines
        props: List of region properties for each spine
    """
    # Extract spine mask
    spine_mask = (segmentation == 1)
    
    # Label connected components
    labeled_spines, num_spines = measure.label(spine_mask, return_num=True, connectivity=1)
    props = measure.regionprops(labeled_spines)
    
    # Filter out small spines
    for prop in props:
        if prop.area < min_size:
            labeled_spines[labeled_spines == prop.label] = 0
    
    # Relabel to make labels consecutive
    labeled_spines, num_spines = measure.label(labeled_spines > 0, return_num=True, connectivity=1)
    props = measure.regionprops(labeled_spines)
    
    return labeled_spines, props

def compute_spine_features(spine_props, frame_idx):
    """
    Compute features for each spine in a timepoint.
    
    Args:
        spine_props: List of region properties
        frame_idx: Index of the current frame
        
    Returns:
        DataFrame with spine features
    """
    features = []
    
    for idx, prop in enumerate(spine_props):
        features.append({
            'frame': frame_idx,
            'label': prop.label,
            'y': prop.centroid[0],
            'x': prop.centroid[1],
            'area': prop.area,
            'bbox': str(prop.bbox),  # Convert to string for CSV compatibility
            'major_axis_length': prop.major_axis_length,
            'minor_axis_length': prop.minor_axis_length,
            'eccentricity': prop.eccentricity,
            'track_id': -1,  # Will be assigned during tracking
            'spine_status': 'unknown'  # Will be assigned during tracking
        })
    
    return pd.DataFrame(features)

def compute_overlap_matrix(labeled_prev, labeled_curr):
    """
    Compute overlap between all spine pairs in consecutive frames.
    
    Args:
        labeled_prev: Labeled spines in previous frame
        labeled_curr: Labeled spines in current frame
        
    Returns:
        overlap_matrix: Matrix of overlap values (range 0-1)
    """
    # Get unique labels (excluding background)
    prev_labels = np.unique(labeled_prev)[1:] if len(np.unique(labeled_prev)) > 1 else []
    curr_labels = np.unique(labeled_curr)[1:] if len(np.unique(labeled_curr)) > 1 else []
    
    if len(prev_labels) == 0 or len(curr_labels) == 0:
        return np.zeros((len(prev_labels), len(curr_labels)))
    
    # Calculate overlap matrix
    overlap_matrix = np.zeros((len(prev_labels), len(curr_labels)))
    
    for i, prev_label in enumerate(prev_labels):
        prev_mask = (labeled_prev == prev_label)
        prev_area = np.sum(prev_mask)
        
        # Dilate previous mask more aggressively to allow for alignment issues
        prev_mask_dilated = prev_mask
        for _ in range(3):
            prev_mask_dilated = binary_dilation(prev_mask_dilated)
        
        for j, curr_label in enumerate(curr_labels):
            curr_mask = (labeled_curr == curr_label)
            curr_area = np.sum(curr_mask)
            
            # Calculate intersection with dilated mask
            intersection = np.sum(np.logical_and(prev_mask_dilated, curr_mask))
            
            # If there's no intersection with the dilated mask, skip
            if intersection == 0:
                continue
            
            # Calculate modified overlap score (more lenient)
            # We'll use a weighted combination of intersection over smallest area
            min_area = min(prev_area, curr_area)
            overlap = 0.7 * (intersection / min_area) + 0.3 * (intersection / (prev_area + curr_area - intersection))
            
            overlap_matrix[i, j] = overlap
    
    return overlap_matrix

def track_spines(frames, min_size=3, max_distance=10, min_overlap=0.3):
    """
    Track spines across all timepoints.
    
    Args:
        frames: 3D array where dimension 0 = timepoints
        min_size: Minimum spine size in pixels
        max_distance: Maximum distance for linking
        min_overlap: Minimum overlap required for linking
        
    Returns:
        all_features: DataFrame with all spine features
        track_labels: Dictionary mapping track ID to track properties
    """
    num_frames = frames.shape[0]
    all_features = []
    next_track_id = 0
    tracks_active = {}  # Currently active tracks
    track_labels = {}   # All tracks with their properties
    
    if num_frames == 0:
        return pd.DataFrame(), {}
    
    # Process first frame
    labeled_spines, props = extract_spine_objects(frames[0], min_size)
    features = compute_spine_features(props, 0)
    
    # Assign new track IDs to all spines in first frame
    for idx, feature in features.iterrows():
        track_id = next_track_id
        features.at[idx, 'track_id'] = track_id
        features.at[idx, 'spine_status'] = 'new'
        tracks_active[track_id] = {
            'last_frame': 0,
            'last_label': feature['label'],
            'positions': [(feature['y'], feature['x'])],
            'areas': [feature['area']],
            'frames': [0]
        }
        track_labels[track_id] = {
            'start_frame': 0,
            'end_frame': 0,
            'status': 'new',
            'duration': 1
        }
        next_track_id += 1
    
    all_features.append(features)
    prev_labeled = labeled_spines
    
    # Process remaining frames
    for frame_idx in range(1, num_frames):
        # Extract spines from current frame
        labeled_spines, props = extract_spine_objects(frames[frame_idx], min_size)
        features = compute_spine_features(props, frame_idx)
        
        if len(props) == 0:
            # No spines in this frame
            all_features.append(features)
            prev_labeled = labeled_spines
            continue
        
        # Compute overlap between previous and current frame
        overlap_matrix = compute_overlap_matrix(prev_labeled, labeled_spines)
        
        # Compute distance between all spine pairs
        prev_centroids = np.array([track['positions'][-1] for track in tracks_active.values()])
        curr_centroids = np.array([prop.centroid for prop in props])
        
        if len(prev_centroids) == 0 or len(curr_centroids) == 0:
            # Handle empty frames
            for j in range(len(props)):
                # Create new track for each spine
                track_id = next_track_id
                next_track_id += 1
                
                spine = props[j]
                tracks_active[track_id] = {
                    'last_frame': frame_idx,
                    'last_label': spine.label,
                    'positions': [(spine.centroid[0], spine.centroid[1])],
                    'areas': [spine.area],
                    'frames': [frame_idx]
                }
                
                track_labels[track_id] = {
                    'start_frame': frame_idx,
                    'end_frame': frame_idx,
                    'status': 'new',
                    'duration': 1
                }
                
                features.at[j, 'track_id'] = track_id
                features.at[j, 'spine_status'] = 'new'
            
            all_features.append(features)
            prev_labeled = labeled_spines
            continue
        
        # Compute distance matrix
        distance_matrix = cdist(prev_centroids, curr_centroids)
        
        # Create cost matrix for assignment
        cost_matrix = np.full((len(tracks_active), len(props)), np.inf)
        
        # Fill in cost matrix based on distance and overlap
        for i, (track_id, track) in enumerate(tracks_active.items()):
            # Find the label index in the previous frame
            last_label = track['last_label']
            last_label_idx = np.where(np.unique(prev_labeled)[1:] == last_label)[0]
            
            for j in range(len(props)):
                # If distance is too large, keep as infinite cost
                distance = distance_matrix[i, j]
                if distance > max_distance:
                    continue
                
                # If we have overlap information, use it
                overlap = overlap_matrix[last_label_idx[0], j] if last_label_idx.size > 0 else 0
                
                # Create a combined cost based on both distance and overlap
                if overlap >= min_overlap or distance <= max_distance * 0.6:
                    # Weight overlap more heavily for close objects, distance for far objects
                    distance_weight = min(1.0, distance / max_distance)
                    overlap_weight = 1.0 - distance_weight
                    
                    # Give more stable spines (ones that have been tracked longer) priority
                    stability_factor = min(1.0, 0.8 + 0.2 * len(track['frames']) / num_frames)
                    
                    # Combined cost: lower is better
                    combined_cost = (distance_weight * (distance / max_distance) + 
                                    overlap_weight * (1.0 - overlap)) / stability_factor
                    
                    cost_matrix[i, j] = combined_cost
            
        # Solve assignment problem
        if np.all(np.isinf(cost_matrix)):
            # All costs are infinite, no valid assignments
            row_indices, col_indices = np.array([]), np.array([])
        else:
            # Replace infinite values with a very large number for the algorithm
            cost_matrix_finite = np.copy(cost_matrix)
            cost_matrix_finite[np.isinf(cost_matrix_finite)] = 1e10
            row_indices, col_indices = linear_sum_assignment(cost_matrix_finite)
            # Filter out assignments with infinite cost
            mask = ~np.isinf(cost_matrix[row_indices, col_indices])
            row_indices, col_indices = row_indices[mask], col_indices[mask]
        
        # Update tracks based on assignments
        assigned_curr = set()
        for i, j in zip(row_indices, col_indices):
            # Get track ID and spine properties
            track_id = list(tracks_active.keys())[i]
            track = tracks_active[track_id]
            spine = props[j]
            
            # Update track
            track['last_frame'] = frame_idx
            track['last_label'] = spine.label
            track['positions'].append((spine.centroid[0], spine.centroid[1]))
            track['areas'].append(spine.area)
            track['frames'].append(frame_idx)
            
            # Update track label
            track_labels[track_id]['end_frame'] = frame_idx
            track_labels[track_id]['duration'] += 1
            track_labels[track_id]['status'] = 'stable'
            
            # Update feature
            features.at[j, 'track_id'] = track_id
            features.at[j, 'spine_status'] = 'existing'
            
            assigned_curr.add(j)
        
        # Handle unassigned current spines (new tracks)
        for j in range(len(props)):
            if j not in assigned_curr:
                # Create new track
                track_id = next_track_id
                next_track_id += 1
                
                spine = props[j]
                tracks_active[track_id] = {
                    'last_frame': frame_idx,
                    'last_label': spine.label,
                    'positions': [(spine.centroid[0], spine.centroid[1])],
                    'areas': [spine.area],
                    'frames': [frame_idx]
                }
                
                track_labels[track_id] = {
                    'start_frame': frame_idx,
                    'end_frame': frame_idx,
                    'status': 'new',
                    'duration': 1
                }
                
                features.at[j, 'track_id'] = track_id
                features.at[j, 'spine_status'] = 'new'
        
        # Handle unassigned previous spines (terminated tracks)
        assigned_prev = set(row_indices)
        for i in range(len(tracks_active)):
            if i not in assigned_prev:
                track_id = list(tracks_active.keys())[i]
                track_labels[track_id]['status'] = 'eliminated'
        
        # Update active tracks with gap closing
        new_tracks_active = {}
        for track_id, track in tracks_active.items():
            # Allow gaps of up to 2 frames for better continuity
            if track['last_frame'] == frame_idx or (frame_idx - track['last_frame']) <= 2:
                # Keep tracks that were active in this frame or just missed up to two frames
                new_tracks_active[track_id] = track
        
        tracks_active = new_tracks_active
        all_features.append(features)
        prev_labeled = labeled_spines
    
    # Combine all features
    all_features_df = pd.concat(all_features, ignore_index=True)
    
    # Update final track status
    for track_id, track in track_labels.items():
        # Calculate stability ratio: what fraction of total frames is this spine present?
        stability_ratio = track['duration'] / num_frames
        
        if stability_ratio >= 0.7:
            # Spine present in 70% or more of frames is considered stable
            track['status'] = 'stable'
        elif track['duration'] == 1:
            if track['start_frame'] == 0:
                track['status'] = 'eliminated'
            elif track['end_frame'] == num_frames - 1:
                track['status'] = 'new'
            else:
                track['status'] = 'transient'
        elif track['duration'] >= 2:
            if track['start_frame'] == 0 and track['end_frame'] == num_frames - 1:
                # Present in first and last frame but not enough frames to be stable
                track['status'] = 'stable' if stability_ratio >= 0.5 else 'transient'
            elif track['start_frame'] > 0 and track['end_frame'] == num_frames - 1:
                track['status'] = 'emerging'
            elif track['start_frame'] == 0 and track['end_frame'] < num_frames - 1:
                track['status'] = 'eliminated'
            else:
                track['status'] = 'transient'
    
    # Update spine status in features dataframe
    for idx, row in all_features_df.iterrows():
        track_id = row['track_id']
        if track_id >= 0:
            status = track_labels[track_id]['status']
            all_features_df.at[idx, 'spine_status'] = status
    
    return all_features_df, track_labels

def create_visualization(segmentations, all_features, track_labels, output_path, metadata=None):
    """
    Create visualization of spine tracking results with day information.
    
    Args:
        segmentations: 3D array of original segmentations
        all_features: DataFrame with all spine features
        track_labels: Dictionary of track properties
        output_path: Path to save visualization
        metadata: Optional metadata dictionary with day information
    """
    num_frames = segmentations.shape[0]
    
    # Create background colormap (black, green, red)
    bg_cmap = ListedColormap(['black', 'green', 'red'])
    
    # Create figure with subplots
    fig, axes = plt.subplots(1, num_frames, figsize=(4*num_frames, 6))
    if num_frames == 1:
        axes = [axes]
    
    # Status colors for legend
    status_colors = {
        'stable': 'blue',
        'new': 'green',
        'eliminated': 'red',
        'emerging': 'cyan',
        'transient': 'yellow'
    }
    
    # Create legend patches
    legend_patches = [mpatches.Patch(color=color, label=status) 
                      for status, color in status_colors.items()]
    
    # Prepare title info with metadata if available
    title_info = ""
    if metadata and metadata.get('animal_id'):
        title_info = f"{metadata.get('animal_id', 'Unknown')}"
        if metadata.get('genotype'):
            title_info = f"{metadata.get('genotype')} {title_info}"
        if metadata.get('segment'):
            title_info = f"{title_info} - {metadata.get('segment')}"
    
    # Plot each frame
    for frame_idx in range(num_frames):
        ax = axes[frame_idx]
        
        # Display background segmentation
        ax.imshow(segmentations[frame_idx], cmap=bg_cmap, alpha=0.3)
        
        # Get spines for this frame
        frame_features = all_features[all_features['frame'] == frame_idx]
        
        # Plot centroids with track colors
        for _, spine in frame_features.iterrows():
            track_id = spine['track_id']
            if track_id >= 0:
                status = spine['spine_status']
                
                # Draw centroid
                ax.plot(spine['x'], spine['y'], 'o', color=status_colors.get(status, 'white'), 
                        markersize=6, alpha=0.8)
                
                # Add track ID label
                ax.text(spine['x'] + 2, spine['y'] + 2, str(track_id), 
                        color='white', fontsize=8, ha='left', va='bottom',
                        bbox=dict(facecolor='black', alpha=0.5, pad=1))
        
        # Set title with day information if available
        if metadata and frame_idx in metadata.get('frame_to_day', {}):
            day_num = metadata['frame_to_day'][frame_idx]
            ax.set_title(f'Day {day_num} (Frame {frame_idx})')
        else:
            ax.set_title(f'Frame {frame_idx}')
        
        ax.set_xlim(0, segmentations.shape[2])
        ax.set_ylim(segmentations.shape[1], 0)  # Invert y-axis
        ax.axis('off')
    
    # Add common legend
    fig.legend(handles=legend_patches, loc='lower center', ncol=len(status_colors),
               bbox_to_anchor=(0.5, 0))
    
    # Add overall title with metadata if available
    if title_info:
        fig.suptitle(title_info, fontsize=14, y=0.98)
        plt.subplots_adjust(top=0.85)  # Make room for title
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)  # Make room for legend
    
    # Save the visualization
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return output_path

def create_visualization_tiff(segmentations, all_features, track_labels, output_path):
    """
    Create a visualization TIFF with color-coded spine tracks.
    
    Args:
        segmentations: 3D array of original segmentations
        all_features: DataFrame with all spine features
        track_labels: Dictionary of track properties
        output_path: Path to save visualization TIFF
    """
    # Create a visualization TIFF with color-coded tracks
    vis_tiff = np.zeros(segmentations.shape, dtype=np.uint8)
    
    # Copy the dendrite (value 2) from original segmentation
    vis_tiff[segmentations == 2] = 2
    
    # For each frame, fill in the spines with their track IDs (shifted by 10 to avoid overlap with dendrite)
    for frame_idx in range(segmentations.shape[0]):
        frame_features = all_features[all_features['frame'] == frame_idx]
        
        for _, spine in frame_features.iterrows():
            track_id = spine['track_id']
            if track_id >= 0:
                # Create binary mask for this spine
                y, x = int(spine['y']), int(spine['x'])
                
                # Get the spine object from the original segmentation
                spine_mask = (segmentations[frame_idx] == 1)
                
                # Find the connected component containing (y, x)
                from skimage import measure
                labeled = measure.label(spine_mask)
                if y < labeled.shape[0] and x < labeled.shape[1] and labeled[y, x] > 0:
                    spine_id = labeled[y, x]
                    spine_mask = (labeled == spine_id)
                    
                    # Assign track ID + 10 to avoid overlap with dendrite (2)
                    vis_tiff[frame_idx][spine_mask] = (track_id % 240) + 10
    
    # Save the visualization TIFF
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tifffile.imwrite(output_path, vis_tiff)
    
    return output_path

def track_spines_with_metadata(segmentation_path, metadata=None, output_dir=None, min_size=3, max_distance=10, min_overlap=0.3):
    """
    Track spines with metadata preservation.
    
    Args:
        segmentation_path: Path to registered segmentation TIFF
        metadata: Metadata dictionary (optional)
        output_dir: Output directory (defaults to the directory of segmentation_path)
        min_size: Minimum spine size in pixels
        max_distance: Maximum distance for linking
        min_overlap: Minimum overlap required for linking
        
    Returns:
        all_features: DataFrame with all spine features
        track_labels: Dictionary mapping track ID to track properties
    """
    # Set up output directory
    if output_dir is None:
        output_dir = os.path.dirname(segmentation_path)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the segmentation
    print(f"Loading {segmentation_path}...")
    segmentations = tifffile.imread(segmentation_path)
    
    # Track spines
    print("Tracking spines...")
    all_features, track_labels = track_spines(
        segmentations,
        min_size=min_size,
        max_distance=max_distance,
        min_overlap=min_overlap
    )
    
    # Calculate tracking statistics
    total_tracks = len(track_labels)
    stable_count = sum(1 for _, track in track_labels.items() if track['status'] == 'stable')
    new_count = sum(1 for _, track in track_labels.items() if track['status'] == 'new' or track['status'] == 'emerging')
    eliminated_count = sum(1 for _, track in track_labels.items() if track['status'] == 'eliminated')
    transient_count = sum(1 for _, track in track_labels.items() if track['status'] == 'transient')
    
    print(f"Tracking complete! Found {total_tracks} total tracks:")
    print(f"  - {stable_count} stable spines")
    print(f"  - {new_count} new/emerging spines")
    print(f"  - {eliminated_count} eliminated spines")
    print(f"  - {transient_count} transient spines")
    
    # Add physical units if available in metadata
    if metadata and 'microns_per_pixel' in metadata:
        microns_per_pixel = metadata['microns_per_pixel']
        all_features['x_um'] = all_features['x'] * microns_per_pixel
        all_features['y_um'] = all_features['y'] * microns_per_pixel
        all_features['area_um2'] = all_features['area'] * (microns_per_pixel ** 2)
    
    # Add day information if available
    if metadata and metadata.get('frame_to_day'):
        all_features['day'] = all_features['frame'].map(
            lambda f: metadata['frame_to_day'].get(f, None)
        )
    
    # Save results
    output_base = os.path.splitext(os.path.basename(segmentation_path))[0]
    csv_path = os.path.join(output_dir, f"{output_base}_spine_tracks.csv")
    all_features.to_csv(csv_path, index=False)
    print(f"Saved spine tracks to {csv_path}")
    
    # Save track summary
    track_summary = pd.DataFrame([
        {
            'track_id': track_id,
            'start_frame': track['start_frame'],
            'end_frame': track['end_frame'],
            'duration': track['duration'],
            'status': track['status'],
        }
        for track_id, track in track_labels.items()
    ])
    
    # Add day information to track summary
    if metadata and metadata.get('frame_to_day'):
        track_summary['start_day'] = track_summary['start_frame'].map(
            lambda f: metadata['frame_to_day'].get(f, None)
        )
        track_summary['end_day'] = track_summary['end_frame'].map(
            lambda f: metadata['frame_to_day'].get(f, None)
        )
    
    summary_path = os.path.join(output_dir, f"{output_base}_track_summary.csv")
    track_summary.to_csv(summary_path, index=False)
    print(f"Saved track summary to {summary_path}")
    
    # Create visualizations
    vis_path = os.path.join(output_dir, f"{output_base}_tracking_visualization.png")
    create_visualization(segmentations, all_features, track_labels, vis_path, metadata=metadata)
    print(f"Saved visualization to {vis_path}")
    
    # Create visualization TIFF
    vis_tiff_path = os.path.join(output_dir, f"{output_base}_tracking_visualization.tif")
    create_visualization_tiff(segmentations, all_features, track_labels, vis_tiff_path)
    print(f"Saved visualization TIFF to {vis_tiff_path}")
    
    # Create a simple tracking report
    report = [
        f"Spine Tracking Report",
        f"===================",
        f"",
        f"Input file: {segmentation_path}",
        f"Number of timepoints: {segmentations.shape[0]}",
        f"",
        f"Tracking parameters:",
        f"  - Minimum spine size: {min_size} pixels",
        f"  - Maximum linking distance: {max_distance} pixels",
        f"  - Minimum overlap: {min_overlap}",
        f"",
        f"Tracking results:",
        f"  - Total tracks: {total_tracks}",
        f"  - Stable spines: {stable_count} ({stable_count/total_tracks*100:.1f}%)",
        f"  - New/emerging spines: {new_count} ({new_count/total_tracks*100:.1f}%)",
        f"  - Eliminated spines: {eliminated_count} ({eliminated_count/total_tracks*100:.1f}%)",
        f"  - Transient spines: {transient_count} ({transient_count/total_tracks*100:.1f}%)",
    ]
    
    # Add day information to report if available
    if metadata and metadata.get('frame_to_day'):
        day_info = [f"  - Frame {frame} = Day {day}" for frame, day in sorted(metadata['frame_to_day'].items())]
        report.extend(["", "Day information:"] + day_info)
    
    report_path = os.path.join(output_dir, f"{output_base}_tracking_report.txt")
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    print(f"Saved tracking report to {report_path}")
    
    return all_features, track_labels