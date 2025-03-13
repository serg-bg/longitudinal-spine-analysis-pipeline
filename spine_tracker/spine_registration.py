#!/usr/bin/env python
"""
Spine Registration Module for Spine Tracking Analysis

This module handles the creation of Maximum Intensity Projections (MIPs) from
3D TIFF stacks and their registration using StackReg. It also provides
functions for assessing alignment quality.
"""

import os
import re
import numpy as np
import tifffile
from skimage import morphology
from skimage.metrics import structural_similarity
from skimage.registration import phase_cross_correlation
import matplotlib.pyplot as plt
import warnings

# Import PyStackReg for registration (matches FIJI's StackReg)
try:
    from pystackreg import StackReg
    HAVE_PYSTACKREG = True
except ImportError:
    HAVE_PYSTACKREG = False
    warnings.warn("PyStackReg not found. This script requires PyStackReg to match FIJI's results. "
                 "Please install it with: 'pip install pystackreg'")

def create_mip(tif_path):
    """
    Create a max-intensity projection from a 3D TIFF stack.
    
    Args:
        tif_path: Path to the TIFF file
        
    Returns:
        MIP image as a 2D numpy array
    """
    # Load the 3D stack
    stack = tifffile.imread(tif_path)
    
    # Check if the stack is 3D (has z dimension)
    if stack.ndim == 3:
        # Create MIP along z-axis (assuming [Z, Y, X] order)
        mip = np.max(stack, axis=0)
    else:
        # If already 2D, just return as is
        mip = stack
    
    # Print diagnostic info
    print(f"    MIP shape: {mip.shape}, dtype: {mip.dtype}, range: [{mip.min()}-{mip.max()}]")
    
    return mip

def pad_to_common_size(images):
    """
    Pad a list of images to a common size (the maximum size found).
    
    Args:
        images: List of 2D numpy arrays
        
    Returns:
        List of padded images
    """
    # Find maximum height and width
    max_height = max(img.shape[0] for img in images)
    max_width = max(img.shape[1] for img in images)
    
    # Pad each image to the maximum size
    padded_images = []
    for img in images:
        h, w = img.shape
        pad_h = max_height - h
        pad_w = max_width - w
        padded_img = np.pad(img, ((0, pad_h), (0, pad_w)), mode='constant')
        padded_images.append(padded_img)
    
    return padded_images

def register_stack_fiji_style(seg_images, num_passes=30):
    """
    Register a stack of segmentation images using StackReg's rigid body transformation,
    closely mimicking FIJI's StackReg workflow with multiple passes.
    
    Args:
        seg_images: List of segmentation images to register
        num_passes: Number of StackReg passes to apply (default: 30)
    
    Returns:
        List of registered segmentation images
    """
    if not HAVE_PYSTACKREG:
        raise ImportError("PyStackReg is required for registration.")
    
    # Convert to uint8 and scale to full range (0-255) - this mimics FIJI's "ScaleConversions" 
    scaled_images = []
    for img in seg_images:
        # Scale the segmentation values: 0->0, 1->127, 2->255
        # This matches what FIJI does with binary images
        scaled = np.zeros_like(img, dtype=np.uint8)
        scaled[img == 1] = 127  # Spines get mid-level value
        scaled[img == 2] = 255  # Dendrite gets max value
        scaled_images.append(scaled)
    
    # Convert to stack for registration
    stack = np.array(scaled_images)
    
    # Initialize StackReg with rigid body transformation (same as FIJI)
    sr = StackReg(StackReg.RIGID_BODY)
    
    # Apply StackReg multiple times in sequence, exactly as done in FIJI
    for pass_num in range(num_passes):
        print(f"    StackReg pass {pass_num+1}/{num_passes}")
        
        # Register the stack using first image as reference
        sr.register_stack(stack, reference='first')
        
        # Transform the stack - no extra parameters
        stack = sr.transform_stack(stack)
    
    # Now convert the aligned stack back to the original segmentation format
    registered_segs = []
    for i in range(len(seg_images)):
        # Create a new segmentation image
        registered = np.zeros_like(seg_images[0])
        
        # Convert back from scaled values to segmentation values
        registered[(stack[i] > 64) & (stack[i] <= 191)] = 1  # Spine
        registered[stack[i] > 191] = 2  # Dendrite
        
        # Clean up interpolation artifacts that create false spine borders
        # This is crucial for segmentation images to maintain discrete objects
        # Method 1: Remove thin/small spine objects that are likely artifacts
        spine_mask = (registered == 1)
        
        # Use morphological opening to remove small/thin artifacts
        cleaned_spine_mask = morphology.opening(spine_mask, morphology.disk(1))
        
        # Remove very small objects (likely interpolation artifacts)
        cleaned_spine_mask = morphology.remove_small_objects(cleaned_spine_mask, min_size=3)
        
        # Apply the cleaned spine mask
        registered_clean = registered.copy()
        registered_clean[spine_mask & ~cleaned_spine_mask] = 0  # Remove artifacts
        
        registered_segs.append(registered_clean)
    
    return registered_segs

def create_multichannel_tif(images, output_path):
    """
    Create a multi-channel TIFF where each channel is a different day.
    Ensures values are preserved exactly for segmentation images.
    
    Args:
        images: List of segmentation images 
        output_path: Path to save the multi-channel TIFF
    """
    # Create all parent directories if they don't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Stack images along a new dimension (channels)
    multichannel = np.stack(images, axis=0).astype(np.uint8)
    
    # Print diagnostic info
    print(f"    Multichannel shape: {multichannel.shape}, range: [{multichannel.min()}-{multichannel.max()}]")
    
    # Save as TIFF with ImageJ-compatible metadata
    tifffile.imwrite(output_path, multichannel, 
                    metadata={'axes': 'CYX'})
    
    print(f"  Created multi-channel TIFF with {len(images)} days as channels: {output_path}")

def assess_alignment_quality(segmentations):
    """
    Assess the alignment quality of each timepoint using the dendrite mask (value 2).
    
    Args:
        segmentations: 3D array of segmentations
        
    Returns:
        quality_scores: List of alignment quality scores for each timepoint
        median_frame_index: Index of the frame with median quality score (good reference)
    """
    num_frames = segmentations.shape[0]
    
    # Extract dendrite masks
    dendrite_masks = [(segmentation == 2).astype(np.uint8) for segmentation in segmentations]
    
    # Calculate pairwise correlation and similarity scores
    similarity_matrix = np.zeros((num_frames, num_frames))
    
    for i in range(num_frames):
        for j in range(i, num_frames):
            if i == j:
                similarity_matrix[i, j] = 1.0
                continue
            
            # Structural similarity between dendrite masks
            ssim_score = structural_similarity(dendrite_masks[i], dendrite_masks[j], 
                                              data_range=1.0)
            
            # Phase cross-correlation to measure shift
            shift, error, _ = phase_cross_correlation(dendrite_masks[i], dendrite_masks[j])
            shift_magnitude = np.sqrt(shift[0]**2 + shift[1]**2)
            
            # Convert shift to a score (smaller shift = higher score)
            max_shift = min(segmentations.shape[1], segmentations.shape[2]) / 5  # 20% of image size
            shift_score = max(0, 1.0 - (shift_magnitude / max_shift))
            
            # Combined score
            combined_score = 0.7 * ssim_score + 0.3 * shift_score
            
            similarity_matrix[i, j] = combined_score
            similarity_matrix[j, i] = combined_score
    
    # Calculate average alignment quality for each timepoint
    quality_scores = np.mean(similarity_matrix, axis=1)
    
    # Identify the median frame (will be used as a reference if needed)
    median_frame_index = np.argsort(quality_scores)[len(quality_scores) // 2]
    
    return quality_scores, median_frame_index

def visualize_alignment_quality(quality_scores, threshold=0.85, metadata=None, output_path=None):
    """
    Create a visualization of alignment quality scores.
    
    Args:
        quality_scores: List of alignment quality scores
        threshold: Quality threshold for good/poor alignment
        metadata: Optional metadata dictionary with day information
        output_path: Path to save visualization
        
    Returns:
        Path to the saved visualization (if output_path is provided)
    """
    plt.figure(figsize=(10, 6))
    bars = plt.bar(range(len(quality_scores)), quality_scores, 
                  color=['green' if score >= threshold else 'red' for score in quality_scores])
    plt.axhline(y=threshold, color='black', linestyle='--', label=f'Threshold ({threshold})')
    
    # Use day information for x-axis labels if available
    if metadata and metadata.get('frame_to_day'):
        x_labels = [f"Day {metadata['frame_to_day'].get(i, '?')}" for i in range(len(quality_scores))]
        plt.xticks(range(len(quality_scores)), x_labels, rotation=45, ha='right')
    
    plt.xlabel('Timepoint')
    plt.ylabel('Alignment Quality Score')
    plt.title('Alignment Quality by Timepoint')
    plt.ylim(0, 1)
    plt.legend()
    
    # Add frame numbers on top of bars
    for i, bar in enumerate(bars):
        # Include day information if available
        if metadata and i in metadata.get('frame_to_day', {}):
            day_info = f"{i}"
        else:
            day_info = f"Frame {i}"
            
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                day_info, ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    # Save the visualization if output path is provided
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved alignment quality visualization to {output_path}")
        return output_path
    
    return None

def process_segment(source_dir, output_dir, num_stackreg_passes=30, quality_threshold=0.85, skip_filter=False):
    """
    Process a segment directory containing 3D TIFF stacks in Segmentation_Labels.
    
    Args:
        source_dir: Path to directory containing Segmentation_Labels
        output_dir: Path to output directory
        num_stackreg_passes: Number of StackReg passes to apply (default: 30)
        quality_threshold: Threshold for alignment quality (default: 0.85)
        skip_filter: Whether to skip alignment quality filtering
        
    Returns:
        Path to registered segmentation TIFF
    """
    # Find Segmentation_Labels directory
    if os.path.basename(source_dir) == 'Segmentation_Labels':
        seg_labels_dir = source_dir
    else:
        seg_labels_dir = os.path.join(source_dir, 'Segmentation_Labels')
        if not os.path.exists(seg_labels_dir):
            raise ValueError(f"Could not find Segmentation_Labels directory in {source_dir}")
    
    # Get all TIFF files in the directory
    tiff_files = [f for f in os.listdir(seg_labels_dir) if f.endswith('.tif') or f.endswith('.tiff')]
    
    # Extract day information from each file using same pattern as MetadataHandler
    day_info = []
    day_pattern = r'Day(\d+)'
    for tiff_file in tiff_files:
        day_match = re.search(day_pattern, tiff_file)
        if day_match:
            day_num = int(day_match.group(1))
            day_info.append((day_num, os.path.join(seg_labels_dir, tiff_file)))
    
    # Sort by day number
    day_info.sort()
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Create MIPs for all timepoints
    print("Creating MIPs for each timepoint...")
    mips = []
    for day, tiff_path in day_info:
        print(f"  Processing Day {day} from {os.path.basename(tiff_path)}")
        mip = create_mip(tiff_path)
        mips.append(mip)
    
    # Step 2: Pad all MIPs to a common size
    print("Padding MIPs to common size...")
    padded_mips = pad_to_common_size(mips)
    
    # Step 3: Register the MIPs
    print(f"Registering MIPs with {num_stackreg_passes} StackReg passes...")
    registered_mips = register_stack_fiji_style(padded_mips, num_passes=num_stackreg_passes)
    
    # Step 4: Create a multichannel TIFF with the registered MIPs
    registered_path = os.path.join(output_dir, 'registered_segmentation.tif')
    create_multichannel_tif(registered_mips, registered_path)
    
    # Step 5: Assess alignment quality
    if not skip_filter:
        print("Assessing alignment quality...")
        # Load the registered segmentation
        segmentations = tifffile.imread(registered_path)
        
        quality_scores, median_frame_index = assess_alignment_quality(segmentations)
        
        # Save alignment quality scores
        quality_data = []
        for i, score in enumerate(quality_scores):
            quality_data.append({
                'timepoint': i,
                'quality_score': score,
                'included': score >= quality_threshold
            })
        
        # Convert to Pandas DataFrame and save as CSV
        import pandas as pd
        quality_df = pd.DataFrame(quality_data)
        quality_csv_path = os.path.join(output_dir, 'alignment_quality.csv')
        quality_df.to_csv(quality_csv_path, index=False)
        print(f"Saved alignment quality scores to {quality_csv_path}")
        
        # Visualize alignment quality
        quality_vis_path = os.path.join(output_dir, 'alignment_quality.png')
        visualize_alignment_quality(quality_scores, threshold=quality_threshold, 
                                  output_path=quality_vis_path)
        
        # Filter out poor quality frames if needed
        good_frames = np.where(quality_scores >= quality_threshold)[0]
        
        if len(good_frames) < len(quality_scores) and len(good_frames) > 0:
            print(f"Filtering out {len(quality_scores) - len(good_frames)} poorly aligned timepoints")
            
            # Extract the good frames
            filtered_segmentations = segmentations[good_frames]
            
            # Save the filtered stack
            filtered_path = os.path.join(output_dir, 'filtered_segmentation.tif')
            tifffile.imwrite(filtered_path, filtered_segmentations)
            print(f"Saved filtered stack to {filtered_path}")
            
            # Return the filtered path instead
            return filtered_path, good_frames
        
    return registered_path, list(range(len(day_info)))