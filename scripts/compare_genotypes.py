#!/usr/bin/env python
"""
Compare Genotypes Script

This script combines data from different genotypes and creates comparative 
visualizations to highlight differences in spine dynamics between genotypes.
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import json
import argparse
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Compare spine dynamics across genotypes.')
    parser.add_argument('--input-dir', '-i', required=True, help='Input directory with processed data')
    parser.add_argument('--output-dir', '-o', required=True, help='Output directory for comparative plots')
    parser.add_argument('--genotypes', '-g', help='Comma-separated list of genotypes to compare (default: all)')
    parser.add_argument('--animals', '-a', help='Comma-separated list of animals to include (default: all)')
    parser.add_argument('--normalize', '-n', action='store_true', help='Normalize spine density across genotypes')
    return parser.parse_args()

def collect_data_by_genotype(input_dir, genotypes=None, animals=None):
    """
    Collect all relevant data organized by genotype.
    
    Args:
        input_dir (str): Path to input directory with processed data
        genotypes (list): Optional list of genotypes to include
        animals (list): Optional list of animals to include
    
    Returns:
        dict: Data organized by genotype and metric
    """
    if genotypes:
        genotype_dirs = [os.path.join(input_dir, g) for g in genotypes]
        genotype_dirs = [d for d in genotype_dirs if os.path.isdir(d)]
    else:
        genotype_dirs = [os.path.join(input_dir, d) for d in os.listdir(input_dir) 
                        if os.path.isdir(os.path.join(input_dir, d))]
    
    print(f"Found {len(genotype_dirs)} genotype directories: {[os.path.basename(d) for d in genotype_dirs]}")
    
    # Initialize data structure
    data = {
        'normalized_spine_density': defaultdict(list),
        'spine_turnover_ratio': defaultdict(list),
        'spine_survival_fraction': defaultdict(list),
        'metadata': defaultdict(list)
    }
    
    # Collect data from each genotype
    for genotype_dir in genotype_dirs:
        genotype = os.path.basename(genotype_dir)
        
        # Get all animal directories
        animal_dirs = [os.path.join(genotype_dir, d) for d in os.listdir(genotype_dir)
                      if os.path.isdir(os.path.join(genotype_dir, d))]
        
        if animals:
            animal_dirs = [d for d in animal_dirs if os.path.basename(d) in animals]
        
        print(f"  Processing genotype {genotype} with {len(animal_dirs)} animals")
        
        # Process each animal
        for animal_dir in animal_dirs:
            animal = os.path.basename(animal_dir)
            
            # Get all segment directories
            segment_dirs = [os.path.join(animal_dir, d) for d in os.listdir(animal_dir)
                           if os.path.isdir(os.path.join(animal_dir, d))]
            
            for segment_dir in segment_dirs:
                segment = os.path.basename(segment_dir)
                
                # Check for validation folder
                validation_dir = os.path.join(segment_dir, 'validation')
                if not os.path.isdir(validation_dir):
                    print(f"    Skipping {animal}/{segment} - no validation directory")
                    continue
                
                # Check for figure_styles folder
                figure_styles_dir = os.path.join(validation_dir, 'figure_styles')
                if not os.path.isdir(figure_styles_dir):
                    print(f"    Skipping {animal}/{segment} - no figure_styles directory")
                    continue
                
                # Load metadata
                metadata_file = os.path.join(segment_dir, 'tracking_metadata.json')
                metadata = None
                if os.path.exists(metadata_file):
                    with open(metadata_file, 'r') as f:
                        metadata = json.load(f)
                        data['metadata'][genotype].append({
                            'animal': animal,
                            'segment': segment,
                            'metadata': metadata
                        })
                
                # Collect normalized spine density data
                density_file = os.path.join(figure_styles_dir, 'normalized_spine_density.csv')
                if os.path.exists(density_file):
                    df = pd.read_csv(density_file)
                    df['animal'] = animal
                    df['segment'] = segment
                    data['normalized_spine_density'][genotype].append(df)
                
                # Collect turnover ratio data
                turnover_file = os.path.join(figure_styles_dir, 'spine_turnover_ratio.csv')
                if os.path.exists(turnover_file):
                    df = pd.read_csv(turnover_file)
                    df['animal'] = animal
                    df['segment'] = segment
                    data['spine_turnover_ratio'][genotype].append(df)
                
                # Collect survival fraction data
                survival_file = os.path.join(figure_styles_dir, 'spine_survival_fraction.csv')
                if os.path.exists(survival_file):
                    df = pd.read_csv(survival_file)
                    df['animal'] = animal
                    df['segment'] = segment
                    data['spine_survival_fraction'][genotype].append(df)
                
                print(f"    Collected data from {animal}/{segment}")
    
    # Combine data frames for each genotype
    for metric in ['normalized_spine_density', 'spine_turnover_ratio', 'spine_survival_fraction']:
        for genotype in data[metric]:
            if data[metric][genotype]:
                data[metric][genotype] = pd.concat(data[metric][genotype], ignore_index=True)
    
    return data

def calculate_actual_spine_density(data, normalize_across_genotypes=False):
    """
    Calculate actual spine density using the raw spine counts.
    Optionally normalize across genotypes.
    
    Args:
        data (dict): Data organized by genotype and metric
        normalize_across_genotypes (bool): Whether to normalize across genotypes
    
    Returns:
        dict: Genotype to DataFrame mapping with actual spine density data
    """
    density_data = {}
    
    # Calculate absolute density for each genotype
    for genotype, df in data['normalized_spine_density'].items():
        if df.empty:
            continue
            
        # Group by day and calculate average spine count
        density_by_day = df.groupby('day')['spine_count'].agg(['mean', 'std', 'count']).reset_index()
        density_by_day['genotype'] = genotype
        density_by_day = density_by_day.rename(columns={'mean': 'spine_count_mean', 
                                                       'std': 'spine_count_std',
                                                       'count': 'segment_count'})
        
        density_data[genotype] = density_by_day
    
    # Combine all genotypes
    all_density = pd.concat(density_data.values(), ignore_index=True)
    
    # Normalize across genotypes if requested
    if normalize_across_genotypes:
        # Get first day values for each genotype
        first_day_values = all_density.groupby('genotype')['spine_count_mean'].first()
        
        # Normalize all values by the first day value for each genotype
        for genotype in density_data:
            first_day_value = first_day_values[genotype]
            density_data[genotype]['spine_density_normalized'] = density_data[genotype]['spine_count_mean'] / first_day_value
    
    return density_data

def calculate_combined_metrics(data):
    """
    Calculate combined metrics for each genotype.
    
    Args:
        data (dict): Data organized by genotype and metric
    
    Returns:
        dict: Combined metrics by genotype
    """
    combined_metrics = {}
    
    for genotype in data['spine_turnover_ratio']:
        turnover_df = data['spine_turnover_ratio'][genotype]
        if turnover_df.empty:
            continue
            
        # Calculate average turnover ratio by day pair
        turnover_by_day = turnover_df.groupby(['day_from', 'day_to'])[['appear', 'disappear', 'persist', 
                                                                     'total_current', 'total_next', 
                                                                     'turnover_ratio']].agg(['mean', 'std', 'count']).reset_index()
        
        # Flatten multi-level columns
        turnover_by_day.columns = ['_'.join(col).strip('_') for col in turnover_by_day.columns.values]
        
        # Add genotype column
        turnover_by_day['genotype'] = genotype
        
        combined_metrics[f"{genotype}_turnover"] = turnover_by_day
    
    for genotype in data['spine_survival_fraction']:
        survival_df = data['spine_survival_fraction'][genotype]
        if survival_df.empty:
            continue
            
        # Calculate average survival fraction by day
        survival_by_day = survival_df.groupby('day')[['surviving_spines', 'survival_fraction']].agg(['mean', 'std', 'count']).reset_index()
        
        # Flatten multi-level columns
        survival_by_day.columns = ['_'.join(col).strip('_') for col in survival_by_day.columns.values]
        
        # Add genotype column
        survival_by_day['genotype'] = genotype
        
        combined_metrics[f"{genotype}_survival"] = survival_by_day
    
    return combined_metrics

def create_comparative_plots(density_data, combined_metrics, output_dir):
    """
    Create comparative plots across genotypes.
    
    Args:
        density_data (dict): Spine density data by genotype
        combined_metrics (dict): Combined metrics by genotype
        output_dir (str): Output directory for plots
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up the style
    sns.set(style="whitegrid")
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.titlesize': 18
    })
    
    # Colors for different genotypes
    colors = {
        'B6': '#1f77b4',  # Blue
        'SRGAP2A_Const': '#ff7f0e',  # Orange
        'Pol2Het': '#2ca02c',  # Green
        'Pol3Het': '#d62728',  # Red
        'Fezf2_Pol2Het': '#9467bd'  # Purple
    }
    
    # Create comparative spine density plot
    if density_data:
        plt.figure(figsize=(12, 8))
        
        # Combine data from all genotypes
        all_density = pd.concat([df.assign(genotype=genotype) for genotype, df in density_data.items()], 
                              ignore_index=True)
        
        # Create plot
        ax = sns.lineplot(data=all_density, x='day', y='spine_count_mean', hue='genotype', 
                       marker='o', linewidth=2, markersize=8, palette=colors)
        
        # Add error bars
        for genotype in all_density['genotype'].unique():
            genotype_data = all_density[all_density['genotype'] == genotype]
            plt.errorbar(genotype_data['day'], genotype_data['spine_count_mean'], 
                       yerr=genotype_data['spine_count_std'], 
                       fmt='none', capsize=5, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        # Customize plot
        plt.title('Spine Density Comparison Across Genotypes')
        plt.xlabel('Day')
        plt.ylabel('Average Spine Count')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(title='Genotype')
        
        # Save plot
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'comparative_spine_density.png'), dpi=300)
        plt.close()
        
        # Save data
        all_density.to_csv(os.path.join(output_dir, 'comparative_spine_density.csv'), index=False)
    
    # Create normalized spine density plot
    if density_data:
        plt.figure(figsize=(12, 8))
        
        # Process data for normalized plot
        normalized_data = []
        for genotype, df in density_data.items():
            # Get first day value
            first_day_value = df.iloc[0]['spine_count_mean'] if not df.empty else 1.0
            
            # Create normalized data
            normalized_df = df.copy()
            normalized_df['normalized_density'] = normalized_df['spine_count_mean'] / first_day_value
            normalized_df['normalized_std'] = normalized_df['spine_count_std'] / first_day_value
            normalized_df['genotype'] = genotype
            
            normalized_data.append(normalized_df)
        
        # Combine data
        all_normalized = pd.concat(normalized_data, ignore_index=True)
        
        # Create plot
        ax = sns.lineplot(data=all_normalized, x='day', y='normalized_density', hue='genotype', 
                       marker='o', linewidth=2, markersize=8, palette=colors)
        
        # Add error bars
        for genotype in all_normalized['genotype'].unique():
            genotype_data = all_normalized[all_normalized['genotype'] == genotype]
            plt.errorbar(genotype_data['day'], genotype_data['normalized_density'], 
                       yerr=genotype_data['normalized_std'], 
                       fmt='none', capsize=5, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        # Customize plot
        plt.title('Normalized Spine Density Comparison Across Genotypes')
        plt.xlabel('Day')
        plt.ylabel('Normalized Spine Density (Day 1 = 1.0)')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(title='Genotype')
        
        # Save plot
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'comparative_normalized_spine_density.png'), dpi=300)
        plt.close()
        
        # Save data
        all_normalized.to_csv(os.path.join(output_dir, 'comparative_normalized_spine_density.csv'), index=False)
    
    # Create turnover ratio comparison plot
    turnover_data = [df for name, df in combined_metrics.items() if 'turnover' in name]
    if turnover_data:
        plt.figure(figsize=(12, 8))
        
        # Combine data
        all_turnover = pd.concat(turnover_data, ignore_index=True)
        
        # Create plot
        ax = sns.lineplot(data=all_turnover, x='day_from', y='turnover_ratio_mean', hue='genotype', 
                       marker='o', linewidth=2, markersize=8, palette=colors)
        
        # Add error bars
        for genotype in all_turnover['genotype'].unique():
            genotype_data = all_turnover[all_turnover['genotype'] == genotype]
            plt.errorbar(genotype_data['day_from'], genotype_data['turnover_ratio_mean'], 
                       yerr=genotype_data['turnover_ratio_std'], 
                       fmt='none', capsize=5, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        # Customize plot
        plt.title('Spine Turnover Ratio Comparison Across Genotypes')
        plt.xlabel('Day')
        plt.ylabel('Average Turnover Ratio')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(title='Genotype')
        
        # Save plot
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'comparative_turnover_ratio.png'), dpi=300)
        plt.close()
        
        # Save data
        all_turnover.to_csv(os.path.join(output_dir, 'comparative_turnover_ratio.csv'), index=False)
    
    # Create survival fraction comparison plot
    survival_data = [df for name, df in combined_metrics.items() if 'survival' in name]
    if survival_data:
        plt.figure(figsize=(12, 8))
        
        # Combine data
        all_survival = pd.concat(survival_data, ignore_index=True)
        
        # Create plot
        ax = sns.lineplot(data=all_survival, x='day', y='survival_fraction_mean', hue='genotype', 
                       marker='o', linewidth=2, markersize=8, palette=colors)
        
        # Add error bars
        for genotype in all_survival['genotype'].unique():
            genotype_data = all_survival[all_survival['genotype'] == genotype]
            plt.errorbar(genotype_data['day'], genotype_data['survival_fraction_mean'], 
                       yerr=genotype_data['survival_fraction_std'], 
                       fmt='none', capsize=5, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        # Customize plot
        plt.title('Spine Survival Fraction Comparison Across Genotypes')
        plt.xlabel('Day')
        plt.ylabel('Average Survival Fraction')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(title='Genotype')
        
        # Save plot
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'comparative_survival_fraction.png'), dpi=300)
        plt.close()
        
        # Save data
        all_survival.to_csv(os.path.join(output_dir, 'comparative_survival_fraction.csv'), index=False)
    
    # Create comprehensive comparison dashboard
    create_dashboard(density_data, combined_metrics, output_dir)

def create_dashboard(density_data, combined_metrics, output_dir):
    """
    Create a comprehensive dashboard with all metrics.
    
    Args:
        density_data (dict): Spine density data by genotype
        combined_metrics (dict): Combined metrics by genotype
        output_dir (str): Output directory for dashboard
    """
    # Set up the style
    sns.set(style="whitegrid")
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 16
    })
    
    # Colors for different genotypes
    colors = {
        'B6': '#1f77b4',  # Blue
        'SRGAP2A_Const': '#ff7f0e',  # Orange
        'Pol2Het': '#2ca02c',  # Green
        'Pol3Het': '#d62728',  # Red
        'Fezf2_Pol2Het': '#9467bd'  # Purple
    }
    
    # Create figure with subplots
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 2, figure=fig, wspace=0.25, hspace=0.3)
    
    # Absolute Spine Density Plot
    ax1 = fig.add_subplot(gs[0, 0])
    if density_data:
        all_density = pd.concat([df.assign(genotype=genotype) for genotype, df in density_data.items()], 
                              ignore_index=True)
        
        sns.lineplot(ax=ax1, data=all_density, x='day', y='spine_count_mean', hue='genotype', 
                   marker='o', linewidth=2, markersize=6, palette=colors)
        
        for genotype in all_density['genotype'].unique():
            genotype_data = all_density[all_density['genotype'] == genotype]
            ax1.errorbar(genotype_data['day'], genotype_data['spine_count_mean'], 
                       yerr=genotype_data['spine_count_std'], 
                       fmt='none', capsize=4, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        ax1.set_title('Absolute Spine Density')
        ax1.set_xlabel('Day')
        ax1.set_ylabel('Average Spine Count')
        ax1.grid(True, linestyle='--', alpha=0.7)
        ax1.legend(title='Genotype')
    
    # Normalized Spine Density Plot
    ax2 = fig.add_subplot(gs[0, 1])
    if density_data:
        normalized_data = []
        for genotype, df in density_data.items():
            first_day_value = df.iloc[0]['spine_count_mean'] if not df.empty else 1.0
            normalized_df = df.copy()
            normalized_df['normalized_density'] = normalized_df['spine_count_mean'] / first_day_value
            normalized_df['normalized_std'] = normalized_df['spine_count_std'] / first_day_value
            normalized_df['genotype'] = genotype
            normalized_data.append(normalized_df)
        
        all_normalized = pd.concat(normalized_data, ignore_index=True)
        
        sns.lineplot(ax=ax2, data=all_normalized, x='day', y='normalized_density', hue='genotype', 
                   marker='o', linewidth=2, markersize=6, palette=colors)
        
        for genotype in all_normalized['genotype'].unique():
            genotype_data = all_normalized[all_normalized['genotype'] == genotype]
            ax2.errorbar(genotype_data['day'], genotype_data['normalized_density'], 
                       yerr=genotype_data['normalized_std'], 
                       fmt='none', capsize=4, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        ax2.set_title('Normalized Spine Density')
        ax2.set_xlabel('Day')
        ax2.set_ylabel('Normalized Density (Day 1 = 1.0)')
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.legend(title='Genotype')
    
    # Turnover Ratio Plot
    ax3 = fig.add_subplot(gs[1, 0])
    turnover_data = [df for name, df in combined_metrics.items() if 'turnover' in name]
    if turnover_data:
        all_turnover = pd.concat(turnover_data, ignore_index=True)
        
        sns.lineplot(ax=ax3, data=all_turnover, x='day_from', y='turnover_ratio_mean', hue='genotype', 
                   marker='o', linewidth=2, markersize=6, palette=colors)
        
        for genotype in all_turnover['genotype'].unique():
            genotype_data = all_turnover[all_turnover['genotype'] == genotype]
            ax3.errorbar(genotype_data['day_from'], genotype_data['turnover_ratio_mean'], 
                       yerr=genotype_data['turnover_ratio_std'], 
                       fmt='none', capsize=4, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        ax3.set_title('Spine Turnover Ratio')
        ax3.set_xlabel('Day')
        ax3.set_ylabel('Average Turnover Ratio')
        ax3.grid(True, linestyle='--', alpha=0.7)
        ax3.legend(title='Genotype')
    
    # Survival Fraction Plot
    ax4 = fig.add_subplot(gs[1, 1])
    survival_data = [df for name, df in combined_metrics.items() if 'survival' in name]
    if survival_data:
        all_survival = pd.concat(survival_data, ignore_index=True)
        
        sns.lineplot(ax=ax4, data=all_survival, x='day', y='survival_fraction_mean', hue='genotype', 
                   marker='o', linewidth=2, markersize=6, palette=colors)
        
        for genotype in all_survival['genotype'].unique():
            genotype_data = all_survival[all_survival['genotype'] == genotype]
            ax4.errorbar(genotype_data['day'], genotype_data['survival_fraction_mean'], 
                       yerr=genotype_data['survival_fraction_std'], 
                       fmt='none', capsize=4, ecolor=colors.get(genotype, 'gray'), alpha=0.5)
        
        ax4.set_title('Spine Survival Fraction')
        ax4.set_xlabel('Day')
        ax4.set_ylabel('Average Survival Fraction')
        ax4.grid(True, linestyle='--', alpha=0.7)
        ax4.legend(title='Genotype')
    
    # Add overall title
    plt.suptitle('Comprehensive Spine Dynamics Comparison Across Genotypes', fontsize=16)
    
    # Save dashboard
    plt.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust for suptitle
    plt.savefig(os.path.join(output_dir, 'spine_dynamics_dashboard.png'), dpi=300)
    plt.close()

def main():
    """Main function."""
    args = parse_arguments()
    
    # Convert comma-separated lists to actual lists
    genotypes = args.genotypes.split(',') if args.genotypes else None
    animals = args.animals.split(',') if args.animals else None
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Collecting data from {args.input_dir}...")
    data = collect_data_by_genotype(args.input_dir, genotypes, animals)
    
    # Check if we have data
    data_available = False
    for metric in data:
        if metric != 'metadata' and any(not df.empty for genotype, df in data[metric].items()):
            data_available = True
            break
    
    if not data_available:
        print("No valid data found. Please check input directory and filters.")
        return
    
    print("Calculating metrics...")
    density_data = calculate_actual_spine_density(data, normalize_across_genotypes=args.normalize)
    combined_metrics = calculate_combined_metrics(data)
    
    print(f"Creating comparative plots in {args.output_dir}...")
    create_comparative_plots(density_data, combined_metrics, args.output_dir)
    
    print("Done!")

if __name__ == "__main__":
    main()