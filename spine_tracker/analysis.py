#!/usr/bin/env python
"""
Analysis Module for Spine Tracking Analysis

This module handles the biological classification of spine tracks and
the generation of literature-style figures with day information.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from skimage import measure
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap

def load_tracking_data(input_dir):
    """
    Load spine tracking data from a directory.
    
    Args:
        input_dir: Directory containing tracking results
        
    Returns:
        spine_tracks: DataFrame with all spine features
        track_summary: DataFrame with track summaries
    """
    # Find spine tracks CSV
    spine_tracks_file = None
    track_summary_file = None
    
    for file in os.listdir(input_dir):
        if file.endswith("_spine_tracks.csv"):
            spine_tracks_file = os.path.join(input_dir, file)
        elif file.endswith("_track_summary.csv"):
            track_summary_file = os.path.join(input_dir, file)
    
    if not spine_tracks_file:
        raise FileNotFoundError(f"No spine tracks CSV found in {input_dir}")
    
    print(f"Loading tracking data from {spine_tracks_file}...")
    spine_tracks = pd.read_csv(spine_tracks_file)
    
    if track_summary_file:
        track_summary = pd.read_csv(track_summary_file)
    else:
        # Create track summary from spine tracks
        track_summary = pd.DataFrame()
        for track_id in spine_tracks['track_id'].unique():
            if track_id < 0:
                continue
                
            track = spine_tracks[spine_tracks['track_id'] == track_id]
            track_summary = track_summary.append({
                'track_id': track_id,
                'start_frame': track['frame'].min(),
                'end_frame': track['frame'].max(),
                'duration': len(track['frame'].unique()),
                'status': track['spine_status'].iloc[0]
            }, ignore_index=True)
    
    return spine_tracks, track_summary

def load_metadata(input_dir):
    """
    Load metadata from a directory.
    
    Args:
        input_dir: Directory containing metadata
        
    Returns:
        metadata: Metadata dictionary
    """
    metadata_file = os.path.join(input_dir, 'tracking_metadata.json')
    
    if not os.path.exists(metadata_file):
        print(f"No metadata file found in {input_dir}")
        return None
    
    from spine_tracker.metadata_handler import MetadataHandler
    metadata_handler = MetadataHandler()
    metadata = metadata_handler.load_metadata(metadata_file)
    
    return metadata

def classify_spine_tracks(spine_tracks, spatial_threshold=1.15, recurrence_window=3):
    """
    Apply biological classification criteria to spine tracks.
    
    Args:
        spine_tracks: DataFrame with spine features
        spatial_threshold: Maximum displacement in microns for spine identity
        recurrence_window: Maximum frames a spine can disappear before being classified as eliminated
        
    Returns:
        enhanced_tracks: DataFrame with biological classification added
    """
    enhanced_tracks = spine_tracks.copy()
    
    # Add physical units if not already present
    if 'x_um' not in enhanced_tracks.columns and 'x' in enhanced_tracks.columns:
        # Default scaling factor (microns per pixel)
        microns_per_pixel = 0.102
        enhanced_tracks['x_um'] = enhanced_tracks['x'] * microns_per_pixel
        enhanced_tracks['y_um'] = enhanced_tracks['y'] * microns_per_pixel
        enhanced_tracks['area_um2'] = enhanced_tracks['area'] * (microns_per_pixel ** 2)
    
    # Add biological status column
    enhanced_tracks['bio_status'] = 'unknown'
    
    # Get all unique track IDs
    track_ids = enhanced_tracks['track_id'].unique()
    track_ids = track_ids[track_ids >= 0]  # Skip unassigned spines
    
    # Analysis tracking spines
    num_frames = enhanced_tracks['frame'].max() + 1
    
    # Check displacements for plausibility
    displacement_errors = []
    spine_displacements = []
    
    for track_id in track_ids:
        track = enhanced_tracks[enhanced_tracks['track_id'] == track_id]
        
        # Calculate displacements in microns
        if 'x_um' in track.columns and len(track) > 1:
            for i in range(len(track) - 1):
                current = track.iloc[i]
                next_point = track.iloc[i+1]
                
                # Calculate displacement in microns
                displacement = np.sqrt((next_point['x_um'] - current['x_um'])**2 + 
                                     (next_point['y_um'] - current['y_um'])**2)
                
                spine_displacements.append({
                    'track_id': track_id,
                    'frame_from': current['frame'],
                    'frame_to': next_point['frame'],
                    'displacement_um': displacement
                })
                
                # Flag large displacements
                if displacement > spatial_threshold:
                    displacement_errors.append({
                        'track_id': track_id,
                        'frame_from': current['frame'],
                        'frame_to': next_point['frame'],
                        'displacement_um': displacement
                    })
        
        # Get frame presence pattern for this track
        frames = set(track['frame'])
        max_frame = max(frames)
        min_frame = min(frames)
        
        # Criterion for stable spines: present in first and last frame
        if 0 in frames and (num_frames - 1) in frames:
            # Check for long gaps
            max_gap = 0
            for i in range(min_frame, max_frame):
                if i not in frames:
                    gap_start = i
                    while i+1 not in frames and i+1 <= max_frame:
                        i += 1
                    gap_end = i
                    gap_length = gap_end - gap_start + 1
                    max_gap = max(max_gap, gap_length)
            
            if max_gap <= recurrence_window:
                # If gaps are small, spine is considered stable
                enhanced_tracks.loc[enhanced_tracks['track_id'] == track_id, 'bio_status'] = 'stable'
            else:
                # If there are long gaps, spine is considered recurrent
                enhanced_tracks.loc[enhanced_tracks['track_id'] == track_id, 'bio_status'] = 'recurrent'
        
        # Criterion for formed (new) spines: not present in first frame but present in last frame
        elif 0 not in frames and (num_frames - 1) in frames:
            enhanced_tracks.loc[enhanced_tracks['track_id'] == track_id, 'bio_status'] = 'formed'
        
        # Criterion for eliminated spines: present in first frame but not in last frame
        elif 0 in frames and (num_frames - 1) not in frames:
            enhanced_tracks.loc[enhanced_tracks['track_id'] == track_id, 'bio_status'] = 'eliminated'
        
        # Other spines are considered transient
        else:
            enhanced_tracks.loc[enhanced_tracks['track_id'] == track_id, 'bio_status'] = 'transient'
    
    # Create DataFrames for tracking errors and displacements
    displacement_errors_df = pd.DataFrame(displacement_errors)
    spine_displacements_df = pd.DataFrame(spine_displacements)
    
    return enhanced_tracks, displacement_errors_df, spine_displacements_df

def analyze_frame_transitions(enhanced_tracks):
    """
    Analyze transitions between frames (spine presence/absence).
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        
    Returns:
        frame_transitions: DataFrame with frame transition statistics
        transition_rates: Dictionary with transition rates
    """
    # Get all unique frames
    frames = sorted(enhanced_tracks['frame'].unique())
    
    # Initialize transition counter
    transition_counts = {
        'appear': [0] * (len(frames) - 1),
        'disappear': [0] * (len(frames) - 1),
        'persist': [0] * (len(frames) - 1),
        'total_spines': [0] * len(frames)
    }
    
    # Count total spines in each frame
    for i, frame in enumerate(frames):
        transition_counts['total_spines'][i] = len(enhanced_tracks[enhanced_tracks['frame'] == frame])
    
    # Count transitions between frames
    for i in range(len(frames) - 1):
        current_frame = frames[i]
        next_frame = frames[i+1]
        
        # Get spines in current and next frame
        current_tracks = set(enhanced_tracks[enhanced_tracks['frame'] == current_frame]['track_id'])
        next_tracks = set(enhanced_tracks[enhanced_tracks['frame'] == next_frame]['track_id'])
        
        # Count spines that persist
        persist = len(current_tracks.intersection(next_tracks))
        transition_counts['persist'][i] = persist
        
        # Count spines that appear
        appear = len(next_tracks - current_tracks)
        transition_counts['appear'][i] = appear
        
        # Count spines that disappear
        disappear = len(current_tracks - next_tracks)
        transition_counts['disappear'][i] = disappear
    
    # Create DataFrame for frame transitions
    frame_transitions = pd.DataFrame({
        'frame_from': frames[:-1],
        'frame_to': frames[1:],
        'total_current': transition_counts['total_spines'][:-1],
        'total_next': transition_counts['total_spines'][1:],
        'persist': transition_counts['persist'],
        'appear': transition_counts['appear'],
        'disappear': transition_counts['disappear']
    })
    
    # Calculate transition rates for each frame transition
    frame_transitions['appear_rate'] = frame_transitions['appear'] / frame_transitions['total_next']
    frame_transitions['disappear_rate'] = frame_transitions['disappear'] / frame_transitions['total_current']
    frame_transitions['persist_rate'] = frame_transitions['persist'] / frame_transitions['total_current']
    
    # Calculate overall rates
    transition_rates = {
        'overall_appear_rate': frame_transitions['appear'].sum() / (frame_transitions['total_next'].sum() + 1e-10),
        'overall_disappear_rate': frame_transitions['disappear'].sum() / (frame_transitions['total_current'].sum() + 1e-10),
        'overall_persist_rate': frame_transitions['persist'].sum() / (frame_transitions['total_current'].sum() + 1e-10),
        'average_appear_rate': frame_transitions['appear_rate'].mean(),
        'average_disappear_rate': frame_transitions['disappear_rate'].mean(),
        'average_persist_rate': frame_transitions['persist_rate'].mean()
    }
    
    return frame_transitions, transition_rates

def create_spine_presence_visualization(enhanced_tracks, output_path, metadata=None):
    """
    Create a visualization of spine presence patterns across frames.
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        output_path: Path to save the visualization
        metadata: Optional metadata dictionary with day information
        
    Returns:
        Path to the saved visualization
    """
    # Get all unique track IDs
    track_ids = enhanced_tracks['track_id'].unique()
    track_ids = track_ids[track_ids >= 0]  # Skip unassigned spines
    
    # Get all unique frames
    frames = sorted(enhanced_tracks['frame'].unique())
    num_frames = len(frames)
    
    # Create a presence matrix: rows = tracks, columns = frames
    presence_matrix = np.zeros((len(track_ids), num_frames))
    
    # Get day labels if available
    if metadata and metadata.get('frame_to_day'):
        from spine_tracker.metadata_handler import MetadataHandler
        metadata_handler = MetadataHandler()
        metadata_handler.metadata = metadata
        day_labels = metadata_handler.get_day_labels()
        frame_labels = [f"Day {metadata['frame_to_day'].get(frame, '?')}" for frame in frames]
    else:
        frame_labels = [f"Frame {frame}" for frame in frames]
    
    # Fill the presence matrix
    for i, track_id in enumerate(track_ids):
        track_frames = set(enhanced_tracks[enhanced_tracks['track_id'] == track_id]['frame'])
        for j, frame in enumerate(frames):
            if frame in track_frames:
                presence_matrix[i, j] = 1
    
    # Sort tracks based on presence pattern for better visualization
    row_sums = presence_matrix.sum(axis=1)
    first_present = np.array([np.argmax(presence_matrix[i, :] > 0) if row_sums[i] > 0 else num_frames 
                             for i in range(len(track_ids))])
    last_present = np.array([num_frames - 1 - np.argmax(presence_matrix[i, ::-1] > 0) if row_sums[i] > 0 else -1 
                            for i in range(len(track_ids))])
    
    # Sort order: first by first frame, then by last frame, then by duration
    sort_order = np.lexsort((row_sums, last_present, first_present))
    
    # Create the plot
    fig = plt.figure(figsize=(12, max(8, len(track_ids) * 0.15)), facecolor='white')
    ax = fig.add_subplot(111)
    ax.imshow(presence_matrix[sort_order, :], aspect='auto', cmap='Blues', 
              interpolation='nearest')
    
    # Set axis labels
    plt.xlabel('Frame')
    plt.ylabel('Spine Track ID')
    
    # Set x-ticks to frame numbers or day labels
    plt.xticks(range(num_frames), frame_labels, rotation=45, ha='right')
    
    # Set y-ticks to sorted track IDs
    sorted_track_ids = track_ids[sort_order]
    if len(sorted_track_ids) <= 30:
        plt.yticks(range(len(sorted_track_ids)), sorted_track_ids)
    else:
        # Too many tracks, show only some tick marks
        tick_indices = np.linspace(0, len(sorted_track_ids) - 1, 30, dtype=int)
        plt.yticks(tick_indices, sorted_track_ids[tick_indices])
    
    # Add title
    plt.title('Spine Presence Patterns Across Frames')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the visualization
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def create_biological_classification_visualization(enhanced_tracks, output_path, metadata=None):
    """
    Create a visualization of biological spine classifications.
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        output_path: Path to save the visualization
        metadata: Optional metadata dictionary with day information
        
    Returns:
        Path to the saved visualization
    """
    # Get all unique frames
    frames = sorted(enhanced_tracks['frame'].unique())
    
    # Count spines by biological status in each frame
    status_counts = []
    for frame in frames:
        frame_tracks = enhanced_tracks[enhanced_tracks['frame'] == frame]
        counts = frame_tracks.groupby('bio_status').size().to_dict()
        counts['frame'] = frame
        status_counts.append(counts)
    
    status_counts_df = pd.DataFrame(status_counts).fillna(0)
    
    # Get day labels if available
    if metadata and metadata.get('frame_to_day'):
        from spine_tracker.metadata_handler import MetadataHandler
        metadata_handler = MetadataHandler()
        metadata_handler.metadata = metadata
        day_labels = metadata_handler.get_day_labels()
        frame_labels = [f"Day {metadata['frame_to_day'].get(frame, '?')}" for frame in frames]
    else:
        frame_labels = [f"Frame {frame}" for frame in frames]
    
    # Create a stacked bar plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Define colors for each status
    colors = {
        'stable': 'blue',
        'formed': 'green',
        'eliminated': 'red',
        'recurrent': 'purple',
        'transient': 'orange'
    }
    
    # Get all possible statuses
    all_statuses = set(enhanced_tracks['bio_status'])
    
    # Create the stacked bar plot
    bottom = np.zeros(len(frames))
    for status in ['stable', 'formed', 'eliminated', 'recurrent', 'transient']:
        if status in status_counts_df.columns:
            values = status_counts_df[status].values
            ax.bar(range(len(frames)), values, bottom=bottom, 
                 label=status.capitalize(), color=colors.get(status, 'gray'))
            bottom += values
    
    # Set axis labels
    ax.set_xlabel('Frame')
    ax.set_ylabel('Number of Spines')
    
    # Set x-ticks to frame numbers or day labels
    ax.set_xticks(range(len(frames)))
    ax.set_xticklabels(frame_labels, rotation=45, ha='right')
    
    # Add legend
    ax.legend()
    
    # Add title
    ax.set_title('Biological Classification of Spines Across Frames')
    
    # Add metadata as text
    if metadata and 'animal_id' in metadata:
        info_text = f"Animal: {metadata.get('animal_id', 'Unknown')}"
        if 'genotype' in metadata:
            info_text = f"Genotype: {metadata['genotype']}, {info_text}"
        if 'segment' in metadata:
            info_text += f", Segment: {metadata['segment']}"
        
        plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=10)
        plt.subplots_adjust(bottom=0.15)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the visualization
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def create_displacement_visualizations(spine_displacements, output_dir, metadata=None):
    """
    Create visualizations for spine displacements.
    
    Args:
        spine_displacements: DataFrame with spine displacement data
        output_dir: Directory to save visualizations
        metadata: Optional metadata dictionary
        
    Returns:
        Dictionary with paths to saved visualizations
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Dictionary to store visualization paths
    vis_paths = {}
    
    # Create displacement histogram (in microns)
    plt.figure(figsize=(10, 6))
    plt.hist(spine_displacements['displacement_um'], bins=20, color='skyblue', edgecolor='black')
    plt.xlabel('Displacement (µm)')
    plt.ylabel('Count')
    plt.title('Spine Displacement Histogram')
    
    # Add vertical line at mean
    mean_displacement = spine_displacements['displacement_um'].mean()
    plt.axvline(mean_displacement, color='red', linestyle='--', 
               label=f'Mean: {mean_displacement:.3f} µm')
    
    # Add metadata as text
    if metadata and 'animal_id' in metadata:
        info_text = f"Animal: {metadata.get('animal_id', 'Unknown')}"
        if 'genotype' in metadata:
            info_text = f"Genotype: {metadata['genotype']}, {info_text}"
        if 'segment' in metadata:
            info_text += f", Segment: {metadata['segment']}"
        
        plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=10)
        plt.subplots_adjust(bottom=0.15)
    
    plt.legend()
    plt.tight_layout()
    
    # Save histogram
    hist_path = os.path.join(output_dir, 'displacement_histogram_um.png')
    plt.savefig(hist_path, dpi=300)
    plt.close()
    vis_paths['histogram'] = hist_path
    
    return vis_paths

def create_normalized_density_plot(enhanced_tracks, output_path, metadata=None):
    """
    Create a normalized spine density plot with day information.
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        output_path: Path to save the plot
        metadata: Optional metadata dictionary with day information
        
    Returns:
        Path to the saved plot
    """
    # Get all unique frames
    frames = sorted(enhanced_tracks['frame'].unique())
    
    # Count spines in each frame
    spine_counts = []
    for frame in frames:
        frame_tracks = enhanced_tracks[enhanced_tracks['frame'] == frame]
        spine_counts.append(len(frame_tracks))
    
    # Normalize by initial spine count
    initial_count = spine_counts[0]
    normalized_density = [count / initial_count for count in spine_counts]
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    
    # Determine x-axis: use days if available, otherwise use sequential frames
    if metadata and metadata.get('frame_to_day'):
        # Use days for x-axis
        days = [metadata['frame_to_day'].get(frame, frame) for frame in frames]
        plt.plot(days, normalized_density, 'o-', color='blue')
        plt.xlabel('Experimental Day')
    else:
        # Use sequential frames as x-axis
        plt.plot(frames, normalized_density, 'o-', color='blue')
        plt.xlabel('Frame')
    
    plt.ylabel('Normalized Spine Density')
    plt.title('Normalized Spine Density Over Time')
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add metadata as text
    if metadata and 'animal_id' in metadata:
        info_text = f"Animal: {metadata.get('animal_id', 'Unknown')}"
        if 'genotype' in metadata:
            info_text = f"Genotype: {metadata['genotype']}, {info_text}"
        if 'segment' in metadata:
            info_text += f", Segment: {metadata['segment']}"
        
        plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=10)
        plt.subplots_adjust(bottom=0.15)
    
    plt.tight_layout()
    
    # Save the plot
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    
    # Create CSV with the data
    csv_path = os.path.splitext(output_path)[0] + '.csv'
    if metadata and metadata.get('frame_to_day'):
        df = pd.DataFrame({
            'frame': frames,
            'day': days,
            'spine_count': spine_counts,
            'normalized_density': normalized_density
        })
    else:
        df = pd.DataFrame({
            'frame': frames,
            'spine_count': spine_counts,
            'normalized_density': normalized_density
        })
    
    df.to_csv(csv_path, index=False)
    
    return output_path

def create_survival_plot(enhanced_tracks, output_path, metadata=None):
    """
    Create a spine survival fraction plot with day information.
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        output_path: Path to save the plot
        metadata: Optional metadata dictionary with day information
        
    Returns:
        Path to the saved plot
    """
    # Get all unique frames
    frames = sorted(enhanced_tracks['frame'].unique())
    
    # Get spines present in the first frame
    initial_tracks = set(enhanced_tracks[enhanced_tracks['frame'] == frames[0]]['track_id'])
    
    # Calculate survival fraction for each subsequent frame
    survival_fraction = []
    for frame in frames:
        frame_tracks = set(enhanced_tracks[enhanced_tracks['frame'] == frame]['track_id'])
        surviving = len(initial_tracks.intersection(frame_tracks))
        survival_fraction.append(surviving / len(initial_tracks))
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    
    # Determine x-axis: use days if available, otherwise use sequential frames
    if metadata and metadata.get('frame_to_day'):
        # Use days for x-axis
        days = [metadata['frame_to_day'].get(frame, frame) for frame in frames]
        plt.plot(days, survival_fraction, 'o-', color='green')
        plt.xlabel('Experimental Day')
    else:
        # Use sequential frames as x-axis
        plt.plot(frames, survival_fraction, 'o-', color='green')
        plt.xlabel('Frame')
    
    plt.ylabel('Survival Fraction')
    plt.title('Spine Survival Fraction Over Time')
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add metadata as text
    if metadata and 'animal_id' in metadata:
        info_text = f"Animal: {metadata.get('animal_id', 'Unknown')}"
        if 'genotype' in metadata:
            info_text = f"Genotype: {metadata['genotype']}, {info_text}"
        if 'segment' in metadata:
            info_text += f", Segment: {metadata['segment']}"
        
        plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=10)
        plt.subplots_adjust(bottom=0.15)
    
    plt.tight_layout()
    
    # Save the plot
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    
    # Create CSV with the data
    csv_path = os.path.splitext(output_path)[0] + '.csv'
    if metadata and metadata.get('frame_to_day'):
        df = pd.DataFrame({
            'frame': frames,
            'day': days,
            'surviving_spines': [s * len(initial_tracks) for s in survival_fraction],
            'survival_fraction': survival_fraction
        })
    else:
        df = pd.DataFrame({
            'frame': frames,
            'surviving_spines': [s * len(initial_tracks) for s in survival_fraction],
            'survival_fraction': survival_fraction
        })
    
    df.to_csv(csv_path, index=False)
    
    return output_path

def create_turnover_plot(enhanced_tracks, frame_transitions, output_path, metadata=None):
    """
    Create a spine turnover ratio plot with day information.
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        frame_transitions: DataFrame with frame transition statistics
        output_path: Path to save the plot
        metadata: Optional metadata dictionary with day information
        
    Returns:
        Path to the saved plot
    """
    # Get all frame transitions
    transitions = frame_transitions.copy()
    
    # Calculate turnover ratio (TOR) for each transition
    # TOR = (formed + eliminated) / 2 * total
    transitions['turnover_ratio'] = (transitions['appear'] + transitions['disappear']) / (2 * transitions['total_current'])
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    
    # Determine x-axis: use days if available, otherwise use sequential frames
    if metadata and metadata.get('frame_to_day'):
        # Use days for x-axis
        days_from = [metadata['frame_to_day'].get(frame, frame) for frame in transitions['frame_from']]
        days_to = [metadata['frame_to_day'].get(frame, frame) for frame in transitions['frame_to']]
        
        # Use average of from/to days as x coordinate
        x_coords = [(days_from[i] + days_to[i]) / 2 for i in range(len(days_from))]
        
        # Create x-tick labels showing the transition
        x_labels = [f"{days_from[i]}-{days_to[i]}" for i in range(len(days_from))]
        
        plt.plot(x_coords, transitions['turnover_ratio'], 'o-', color='red')
        plt.xticks(x_coords, x_labels, rotation=45, ha='right')
        plt.xlabel('Days Transition')
    else:
        # Use midpoint of frame transitions as x coordinate
        x_coords = [(transitions['frame_from'] + transitions['frame_to']) / 2]
        
        # Create x-tick labels showing the transition
        x_labels = [f"{int(frame_from)}-{int(frame_to)}" for frame_from, frame_to in zip(transitions['frame_from'], transitions['frame_to'])]
        
        plt.plot(x_coords, transitions['turnover_ratio'], 'o-', color='red')
        plt.xticks(x_coords, x_labels, rotation=45, ha='right')
        plt.xlabel('Frame Transition')
    
    plt.ylabel('Turnover Ratio')
    plt.title('Spine Turnover Ratio Between Timepoints')
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add metadata as text
    if metadata and 'animal_id' in metadata:
        info_text = f"Animal: {metadata.get('animal_id', 'Unknown')}"
        if 'genotype' in metadata:
            info_text = f"Genotype: {metadata['genotype']}, {info_text}"
        if 'segment' in metadata:
            info_text += f", Segment: {metadata['segment']}"
        
        plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=10)
        plt.subplots_adjust(bottom=0.15)
    
    plt.tight_layout()
    
    # Save the plot
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    
    # Create CSV with the data
    csv_path = os.path.splitext(output_path)[0] + '.csv'
    if metadata and metadata.get('frame_to_day'):
        transitions['day_from'] = days_from
        transitions['day_to'] = days_to
    
    # Select relevant columns for CSV
    csv_columns = ['frame_from', 'frame_to', 'appear', 'disappear', 'persist', 
                  'total_current', 'total_next', 'turnover_ratio']
    
    if 'day_from' in transitions.columns:
        csv_columns = ['frame_from', 'frame_to', 'day_from', 'day_to'] + csv_columns[2:]
    
    transitions[csv_columns].to_csv(csv_path, index=False)
    
    return output_path

def create_literature_style_figure(enhanced_tracks, frame_transitions, output_path, metadata=None):
    """
    Create a combined literature-style figure with all key metrics.
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        frame_transitions: DataFrame with frame transition statistics
        output_path: Path to save the figure
        metadata: Optional metadata dictionary with day information
        
    Returns:
        Path to the saved figure
    """
    # Create a figure with 3 subplots (normalized density, survival, turnover)
    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
    
    # Get all unique frames
    frames = sorted(enhanced_tracks['frame'].unique())
    
    # Determine x-axis: use days if available, otherwise use sequential frames
    if metadata and metadata.get('frame_to_day'):
        # Use days for x-axis
        days = [metadata['frame_to_day'].get(frame, frame) for frame in frames]
        x_axis = days
        x_label = 'Experimental Day'
        
        # For transitions
        days_from = [metadata['frame_to_day'].get(frame, frame) for frame in frame_transitions['frame_from']]
        days_to = [metadata['frame_to_day'].get(frame, frame) for frame in frame_transitions['frame_to']]
        transition_x = [(days_from[i] + days_to[i]) / 2 for i in range(len(days_from))]
        transition_labels = [f"{days_from[i]}-{days_to[i]}" for i in range(len(days_from))]
    else:
        # Use sequential frames
        x_axis = frames
        x_label = 'Frame'
        
        # For transitions
        transition_x = [(frame_transitions['frame_from'] + frame_transitions['frame_to']) / 2]
        transition_labels = [f"{int(frame_from)}-{int(frame_to)}" for frame_from, frame_to in zip(frame_transitions['frame_from'], frame_transitions['frame_to'])]
    
    # ------------------------------
    # 1. Normalized Spine Density
    # ------------------------------
    ax1 = plt.subplot(gs[0])
    
    # Count spines in each frame
    spine_counts = []
    for frame in frames:
        frame_tracks = enhanced_tracks[enhanced_tracks['frame'] == frame]
        spine_counts.append(len(frame_tracks))
    
    # Normalize by initial spine count
    initial_count = spine_counts[0]
    normalized_density = [count / initial_count for count in spine_counts]
    
    ax1.plot(x_axis, normalized_density, 'o-', color='blue', linewidth=2)
    ax1.set_ylabel('Normalized Spine Density')
    ax1.set_title('A. Normalized Spine Density Over Time')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # ------------------------------
    # 2. Survival Fraction
    # ------------------------------
    ax2 = plt.subplot(gs[1])
    
    # Get spines present in the first frame
    initial_tracks = set(enhanced_tracks[enhanced_tracks['frame'] == frames[0]]['track_id'])
    
    # Calculate survival fraction for each frame
    survival_fraction = []
    for frame in frames:
        frame_tracks = set(enhanced_tracks[enhanced_tracks['frame'] == frame]['track_id'])
        surviving = len(initial_tracks.intersection(frame_tracks))
        survival_fraction.append(surviving / len(initial_tracks))
    
    ax2.plot(x_axis, survival_fraction, 'o-', color='green', linewidth=2)
    ax2.set_ylabel('Survival Fraction')
    ax2.set_title('B. Spine Survival Fraction Over Time')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # ------------------------------
    # 3. Turnover Ratio
    # ------------------------------
    ax3 = plt.subplot(gs[2])
    
    # Calculate turnover ratio
    frame_transitions['turnover_ratio'] = (frame_transitions['appear'] + frame_transitions['disappear']) / (2 * frame_transitions['total_current'])
    
    ax3.plot(transition_x, frame_transitions['turnover_ratio'], 'o-', color='red', linewidth=2)
    ax3.set_ylabel('Turnover Ratio')
    ax3.set_title('C. Spine Turnover Ratio Between Timepoints')
    ax3.set_xticks(transition_x)
    ax3.set_xticklabels(transition_labels, rotation=45, ha='right')
    ax3.grid(True, linestyle='--', alpha=0.7)
    
    # ------------------------------
    # Common Formatting
    # ------------------------------
    for ax in [ax1, ax2]:
        ax.set_xticks(x_axis)
        ax.set_xticklabels(x_axis)
    
    ax3.set_xlabel(x_label)
    
    # Add metadata as text
    if metadata and 'animal_id' in metadata:
        info_text = f"Animal: {metadata.get('animal_id', 'Unknown')}"
        if 'genotype' in metadata:
            info_text = f"Genotype: {metadata['genotype']}, {info_text}"
        if 'segment' in metadata:
            info_text += f", Segment: {metadata['segment']}"
        
        plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=12, weight='bold')
        plt.subplots_adjust(bottom=0.1)
    
    plt.tight_layout()
    
    # Save the figure
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create CSV with the combined data
    csv_path = os.path.splitext(output_path)[0] + '.csv'
    
    # Create a dictionary with all data
    combined_data = {
        'frame': frames,
        'spine_count': spine_counts,
        'normalized_density': normalized_density,
        'survival_fraction': survival_fraction
    }
    
    if metadata and metadata.get('frame_to_day'):
        combined_data['day'] = days
    
    # Add turnover data aligned with frames
    turnover_data = {}
    for i in range(len(frames) - 1):
        current_frame = frames[i]
        next_frame = frames[i + 1]
        
        # Find the corresponding row in frame_transitions
        transition_row = frame_transitions[(frame_transitions['frame_from'] == current_frame) & 
                                         (frame_transitions['frame_to'] == next_frame)]
        
        if not transition_row.empty:
            turnover_data[current_frame] = transition_row['turnover_ratio'].values[0]
    
    # Add turnover data to combined data where available
    combined_data['turnover_ratio'] = [turnover_data.get(frame, np.nan) for frame in frames]
    
    # Create DataFrame and save to CSV
    combined_df = pd.DataFrame(combined_data)
    combined_df.to_csv(csv_path, index=False)
    
    return output_path

def create_frame_transition_report(frame_transitions, output_path):
    """
    Create a report of frame transitions with statistics.
    
    Args:
        frame_transitions: DataFrame with frame transition statistics
        output_path: Path to save the report
        
    Returns:
        Path to the saved report
    """
    # Calculate statistics for the report
    stats = {
        'avg_appear_rate': frame_transitions['appear_rate'].mean(),
        'std_appear_rate': frame_transitions['appear_rate'].std(),
        'avg_disappear_rate': frame_transitions['disappear_rate'].mean(),
        'std_disappear_rate': frame_transitions['disappear_rate'].std(),
        'avg_persist_rate': frame_transitions['persist_rate'].mean(),
        'std_persist_rate': frame_transitions['persist_rate'].std(),
        'avg_turnover_ratio': frame_transitions['turnover_ratio'].mean() if 'turnover_ratio' in frame_transitions.columns else np.nan,
        'std_turnover_ratio': frame_transitions['turnover_ratio'].std() if 'turnover_ratio' in frame_transitions.columns else np.nan
    }
    
    # Create a DataFrame for the statistics
    stats_df = pd.DataFrame({
        'metric': ['Appearance Rate', 'Disappearance Rate', 'Persistence Rate', 'Turnover Ratio'],
        'mean': [stats['avg_appear_rate'], stats['avg_disappear_rate'], stats['avg_persist_rate'], stats['avg_turnover_ratio']],
        'std': [stats['std_appear_rate'], stats['std_disappear_rate'], stats['std_persist_rate'], stats['std_turnover_ratio']]
    })
    
    # Save the statistics to CSV
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    stats_df.to_csv(output_path, index=False)
    
    return output_path

def create_frame_transition_visualization(frame_transitions, output_path, metadata=None):
    """
    Create a visualization of frame transition rates.
    
    Args:
        frame_transitions: DataFrame with frame transition statistics
        output_path: Path to save the visualization
        metadata: Optional metadata dictionary with day information
        
    Returns:
        Path to the saved visualization
    """
    # Create figure with 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    
    # Determine x-axis: use days if available, otherwise use sequential frames
    if metadata and metadata.get('frame_to_day'):
        # Use days for x-axis
        days_from = [metadata['frame_to_day'].get(frame, frame) for frame in frame_transitions['frame_from']]
        days_to = [metadata['frame_to_day'].get(frame, frame) for frame in frame_transitions['frame_to']]
        
        # Use average of from/to days as x coordinate
        x_coords = [(days_from[i] + days_to[i]) / 2 for i in range(len(days_from))]
        
        # Create x-tick labels showing the transition
        x_labels = [f"{days_from[i]}-{days_to[i]}" for i in range(len(days_from))]
        
        x_label = 'Days Transition'
    else:
        # Use midpoint of frame transitions as x coordinate
        x_coords = [(frame_transitions['frame_from'] + frame_transitions['frame_to']) / 2]
        
        # Create x-tick labels showing the transition
        x_labels = [f"{int(frame_from)}-{int(frame_to)}" for frame_from, frame_to in zip(frame_transitions['frame_from'], frame_transitions['frame_to'])]
        
        x_label = 'Frame Transition'
    
    # Plot appearance rate
    ax1.plot(x_coords, frame_transitions['appear_rate'], 'o-', color='green')
    ax1.set_ylabel('Appearance Rate')
    ax1.set_title('A. Spine Appearance Rate')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_xticks(x_coords)
    ax1.set_xticklabels(x_labels, rotation=45, ha='right')
    
    # Plot disappearance rate
    ax2.plot(x_coords, frame_transitions['disappear_rate'], 'o-', color='red')
    ax2.set_ylabel('Disappearance Rate')
    ax2.set_title('B. Spine Disappearance Rate')
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_xticks(x_coords)
    ax2.set_xticklabels(x_labels, rotation=45, ha='right')
    
    # Plot persistence rate
    ax3.plot(x_coords, frame_transitions['persist_rate'], 'o-', color='blue')
    ax3.set_ylabel('Persistence Rate')
    ax3.set_title('C. Spine Persistence Rate')
    ax3.grid(True, linestyle='--', alpha=0.7)
    ax3.set_xticks(x_coords)
    ax3.set_xticklabels(x_labels, rotation=45, ha='right')
    ax3.set_xlabel(x_label)
    
    # Add metadata as text
    if metadata and 'animal_id' in metadata:
        info_text = f"Animal: {metadata.get('animal_id', 'Unknown')}"
        if 'genotype' in metadata:
            info_text = f"Genotype: {metadata['genotype']}, {info_text}"
        if 'segment' in metadata:
            info_text += f", Segment: {metadata['segment']}"
        
        plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=10)
        plt.subplots_adjust(bottom=0.1)
    
    plt.tight_layout()
    
    # Save the visualization
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def create_biological_summary(enhanced_tracks, output_path, metadata=None):
    """
    Create a summary of biological spine classifications.
    
    Args:
        enhanced_tracks: DataFrame with enhanced spine tracks
        output_path: Path to save the summary
        metadata: Optional metadata dictionary
        
    Returns:
        Path to the saved summary
    """
    # Count spines by biological status
    status_counts = enhanced_tracks.groupby('bio_status')['track_id'].nunique()
    
    # Calculate percentages
    total_tracks = status_counts.sum()
    status_percentages = (status_counts / total_tracks * 100).round(1)
    
    # Create summary DataFrame
    summary = pd.DataFrame({
        'status': status_counts.index,
        'count': status_counts.values,
        'percentage': status_percentages.values
    })
    
    # Add experiment information if available
    if metadata:
        summary['animal_id'] = metadata.get('animal_id', 'Unknown')
        summary['genotype'] = metadata.get('genotype', 'Unknown')
        summary['segment'] = metadata.get('segment', 'Unknown')
    
    # Save the summary
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    summary.to_csv(output_path, index=False)
    
    return output_path

def create_html_report(input_dir, output_path, enhanced_tracks=None, metadata=None):
    """
    Create an HTML report with all analysis results.
    
    Args:
        input_dir: Directory containing analysis results
        output_path: Path to save the HTML report
        enhanced_tracks: DataFrame with enhanced spine tracks (optional)
        metadata: Metadata dictionary (optional)
        
    Returns:
        Path to the saved HTML report
    """
    # Load data if not provided
    if enhanced_tracks is None:
        spine_tracks, track_summary = load_tracking_data(input_dir)
        
        # Try to find enhanced tracks
        enhanced_tracks_path = os.path.join(input_dir, 'enhanced_spine_tracks.csv')
        if os.path.exists(enhanced_tracks_path):
            enhanced_tracks = pd.read_csv(enhanced_tracks_path)
        else:
            enhanced_tracks = spine_tracks
    
    if metadata is None:
        metadata = load_metadata(input_dir)
    
    # Create HTML report
    html_content = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "<title>Spine Tracking Analysis Report</title>",
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 20px; }",
        "h1, h2, h3 { color: #2c3e50; }",
        "img { max-width: 100%; height: auto; border: 1px solid #ddd; margin: 10px 0; }",
        ".summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }",
        "table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }",
        "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
        "th { background-color: #f2f2f2; }",
        "tr:nth-child(even) { background-color: #f9f9f9; }",
        "</style>",
        "</head>",
        "<body>"
    ]
    
    # Add title
    html_content.append("<h1>Spine Tracking Analysis Report</h1>")
    
    # Add metadata section
    html_content.append("<h2>Experiment Information</h2>")
    html_content.append("<div class='summary'>")
    
    if metadata:
        html_content.append("<table>")
        html_content.append("<tr><th>Parameter</th><th>Value</th></tr>")
        
        if 'animal_id' in metadata:
            html_content.append(f"<tr><td>Animal ID</td><td>{metadata['animal_id']}</td></tr>")
        
        if 'genotype' in metadata:
            html_content.append(f"<tr><td>Genotype</td><td>{metadata['genotype']}</td></tr>")
        
        if 'segment' in metadata:
            html_content.append(f"<tr><td>Segment</td><td>{metadata['segment']}</td></tr>")
        
        if 'timepoints' in metadata:
            html_content.append(f"<tr><td>Number of Timepoints</td><td>{len(metadata['timepoints'])}</td></tr>")
        
        if 'removed_days' in metadata and metadata['removed_days']:
            html_content.append(f"<tr><td>Removed Days</td><td>{', '.join(map(str, metadata['removed_days']))}</td></tr>")
        
        html_content.append("</table>")
        
        # Add day mapping if available
        if 'frame_to_day' in metadata and metadata['frame_to_day']:
            html_content.append("<h3>Frame to Day Mapping</h3>")
            html_content.append("<table>")
            html_content.append("<tr><th>Frame</th><th>Experimental Day</th></tr>")
            
            for frame, day in sorted(metadata['frame_to_day'].items()):
                html_content.append(f"<tr><td>{frame}</td><td>{day}</td></tr>")
            
            html_content.append("</table>")
    else:
        html_content.append("<p>No metadata available.</p>")
    
    html_content.append("</div>")
    
    # Add spine tracking summary
    html_content.append("<h2>Spine Tracking Summary</h2>")
    html_content.append("<div class='summary'>")
    
    # Count spines by biological status
    if 'bio_status' in enhanced_tracks.columns:
        status_counts = enhanced_tracks.groupby('bio_status')['track_id'].nunique()
        total_tracks = status_counts.sum()
        
        html_content.append("<table>")
        html_content.append("<tr><th>Spine Type</th><th>Count</th><th>Percentage</th></tr>")
        
        for status, count in status_counts.items():
            percentage = count / total_tracks * 100
            html_content.append(f"<tr><td>{status.capitalize()}</td><td>{count}</td><td>{percentage:.1f}%</td></tr>")
        
        html_content.append(f"<tr><td><strong>Total</strong></td><td>{total_tracks}</td><td>100.0%</td></tr>")
        html_content.append("</table>")
    else:
        html_content.append("<p>No biological classification data available.</p>")
    
    html_content.append("</div>")
    
    # Add visualizations section
    html_content.append("<h2>Visualizations</h2>")
    
    # List all visualization files in the input directory
    vis_files = {
        'density': 'normalized_spine_density.png',
        'survival': 'spine_survival_fraction.png',
        'turnover': 'spine_turnover_ratio.png',
        'combined': 'combined_literature_style_figure.png',
        'bio_class': 'biological_classification.png',
        'presence': 'spine_presence_patterns.png',
        'transitions': 'frame_transition_rates.png'
    }
    
    # Add each visualization if it exists
    for vis_type, vis_file in vis_files.items():
        vis_path = os.path.join(input_dir, vis_file)
        if os.path.exists(vis_path):
            # Determine appropriate title
            if vis_type == 'density':
                title = "Normalized Spine Density"
            elif vis_type == 'survival':
                title = "Spine Survival Fraction"
            elif vis_type == 'turnover':
                title = "Spine Turnover Ratio"
            elif vis_type == 'combined':
                title = "Combined Literature-Style Figure"
            elif vis_type == 'bio_class':
                title = "Biological Classification"
            elif vis_type == 'presence':
                title = "Spine Presence Patterns"
            elif vis_type == 'transitions':
                title = "Frame Transition Rates"
            else:
                title = vis_file.replace('_', ' ').replace('.png', '').title()
            
            html_content.append(f"<h3>{title}</h3>")
            
            # Create a relative path for the image
            rel_path = os.path.relpath(vis_path, os.path.dirname(output_path))
            html_content.append(f"<img src='{rel_path}' alt='{title}'>")
    
    # Close HTML
    html_content.append("</body>")
    html_content.append("</html>")
    
    # Save the HTML report
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('\n'.join(html_content))
    
    return output_path

def process_tracking_results(input_dir, output_dir=None, spatial_threshold=1.15, recurrence_window=3):
    """
    Process tracking results to generate biological classifications and visualizations.
    
    Args:
        input_dir: Directory containing tracking results
        output_dir: Output directory (defaults to input_dir)
        spatial_threshold: Maximum displacement in microns for spine identity
        recurrence_window: Maximum frames a spine can disappear before being classified as eliminated
        
    Returns:
        Dictionary with paths to generated files
    """
    # Set up output directory
    if output_dir is None:
        output_dir = input_dir
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load tracking data and metadata
    spine_tracks, track_summary = load_tracking_data(input_dir)
    metadata = load_metadata(input_dir)
    
    print(f"Loaded {len(spine_tracks)} spine detections across {len(spine_tracks['track_id'].unique()) - 1} tracks")
    
    # Apply biological classification
    print("Applying biological criteria...")
    enhanced_tracks, displacement_errors, spine_displacements = classify_spine_tracks(
        spine_tracks, 
        spatial_threshold=spatial_threshold, 
        recurrence_window=recurrence_window
    )
    
    # Save enhanced tracks
    enhanced_tracks_path = os.path.join(output_dir, 'enhanced_spine_tracks.csv')
    enhanced_tracks.to_csv(enhanced_tracks_path, index=False)
    print(f"Saved enhanced tracks to {enhanced_tracks_path}")
    
    # Create output subdirectories
    bio_class_dir = os.path.join(output_dir, 'biological_classification')
    os.makedirs(bio_class_dir, exist_ok=True)
    
    reg_metrics_dir = os.path.join(output_dir, 'registration_metrics')
    os.makedirs(reg_metrics_dir, exist_ok=True)
    
    figure_dir = os.path.join(output_dir, 'figure_styles')
    os.makedirs(figure_dir, exist_ok=True)
    
    # Dictionary to store all output paths
    outputs = {
        'enhanced_tracks': enhanced_tracks_path
    }
    
    # Save displacement data
    if len(displacement_errors) > 0:
        print(f"Found {len(displacement_errors)} tracks with biologically implausible displacements")
        errors_path = os.path.join(reg_metrics_dir, 'displacement_errors.csv')
        pd.DataFrame(displacement_errors).to_csv(errors_path, index=False)
        outputs['displacement_errors'] = errors_path
    
    displacements_path = os.path.join(reg_metrics_dir, 'spine_displacements.csv')
    pd.DataFrame(spine_displacements).to_csv(displacements_path, index=False)
    outputs['spine_displacements'] = displacements_path
    
    # Calculate displacement statistics in microns
    if 'displacement_um' in pd.DataFrame(spine_displacements).columns:
        displacements_um = pd.DataFrame(spine_displacements)['displacement_um']
        displacements_um_path = os.path.join(reg_metrics_dir, 'spine_displacements_um.csv')
        pd.DataFrame({'displacement_um': displacements_um}).to_csv(displacements_um_path, index=False)
        outputs['spine_displacements_um'] = displacements_um_path
        
        # Print summary statistics
        print("Displacement Statistics (microns):")
        print(f"Mean: {displacements_um.mean():.3f} μm")
        print(f"Median: {displacements_um.median():.3f} μm")
        print(f"Max: {displacements_um.max():.3f} μm")
        print(f"Std: {displacements_um.std():.3f} μm")
        
        # Create displacement histogram
        hist_path = os.path.join(reg_metrics_dir, 'displacement_histogram_um.png')
        outputs.update(create_displacement_visualizations(pd.DataFrame(spine_displacements), reg_metrics_dir, metadata))
    
    # Create biological classification visualization
    print("Creating biological classification visualizations...")
    bio_class_path = os.path.join(bio_class_dir, 'biological_classification.png')
    outputs['bio_classification'] = create_biological_classification_visualization(enhanced_tracks, bio_class_path, metadata)
    
    # Create spine presence patterns visualization
    presence_path = os.path.join(bio_class_dir, 'spine_presence_patterns.png')
    outputs['presence_patterns'] = create_spine_presence_visualization(enhanced_tracks, presence_path, metadata)
    
    # Analyze frame transitions
    print("Analyzing frame-to-frame transitions...")
    frame_transitions, transition_rates = analyze_frame_transitions(enhanced_tracks)
    
    # Save frame transition data
    transitions_path = os.path.join(bio_class_dir, 'frame_transitions.csv')
    frame_transitions.to_csv(transitions_path, index=False)
    outputs['frame_transitions'] = transitions_path
    
    # Create frame transition visualization
    transition_vis_path = os.path.join(bio_class_dir, 'frame_transition_rates.png')
    outputs['transition_rates_vis'] = create_frame_transition_visualization(frame_transitions, transition_vis_path, metadata)
    
    # Save transition statistics
    transition_stats_path = os.path.join(bio_class_dir, 'frame_transition_statistics.csv')
    outputs['transition_stats'] = create_frame_transition_report(frame_transitions, transition_stats_path)
    
    # Create biological summary
    bio_summary_path = os.path.join(bio_class_dir, 'biological_summary.csv')
    outputs['bio_summary'] = create_biological_summary(enhanced_tracks, bio_summary_path, metadata)
    
    # Create literature-style figures
    print("Creating figure style plots...")
    
    # Normalized spine density
    density_path = os.path.join(figure_dir, 'normalized_spine_density.png')
    outputs['density_plot'] = create_normalized_density_plot(enhanced_tracks, density_path, metadata)
    
    # Spine survival fraction
    survival_path = os.path.join(figure_dir, 'spine_survival_fraction.png')
    outputs['survival_plot'] = create_survival_plot(enhanced_tracks, survival_path, metadata)
    
    # Spine turnover ratio
    turnover_path = os.path.join(figure_dir, 'spine_turnover_ratio.png')
    outputs['turnover_plot'] = create_turnover_plot(enhanced_tracks, frame_transitions, turnover_path, metadata)
    
    # Combined literature-style figure
    combined_path = os.path.join(figure_dir, 'combined_literature_style_figure.png')
    outputs['combined_figure'] = create_literature_style_figure(enhanced_tracks, frame_transitions, combined_path, metadata)
    
    # Create HTML report
    print("Creating HTML report...")
    report_path = os.path.join(output_dir, 'biological_criteria_report.html')
    outputs['html_report'] = create_html_report(output_dir, report_path, enhanced_tracks, metadata)
    
    return outputs