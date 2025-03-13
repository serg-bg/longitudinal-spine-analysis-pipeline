#!/usr/bin/env python
"""
Main entry point for the Spine Tracking Analysis Pipeline.

This script provides a command-line interface to the various pipeline components.
"""

import os
import sys
import argparse
import subprocess

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Spine Tracking Analysis Pipeline')
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Test command
    test_parser = subparsers.add_parser('test', help='Run the test workflow on a single segment')
    
    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Process multiple segments in batch')
    batch_parser.add_argument('--genotypes', '-g', help='Comma-separated list of genotypes to process')
    batch_parser.add_argument('--animals', '-a', help='Comma-separated list of animals to process')
    
    # Dry run command
    dry_run_parser = subparsers.add_parser('dry-run', help='Count segments without processing')
    dry_run_parser.add_argument('--genotypes', '-g', help='Comma-separated list of genotypes to process')
    
    return parser.parse_args()

def run_script(script_name, **kwargs):
    """Run a shell script with given arguments."""
    script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts')
    script_path = os.path.join(script_dir, script_name)
    
    # Build command with arguments
    command = [script_path]
    for key, value in kwargs.items():
        if value is not None:
            command.extend([f'--{key.replace("_", "-")}', value])
    
    # Run the script
    process = subprocess.Popen(command)
    process.wait()
    
    return process.returncode

def main():
    """Main function to dispatch commands."""
    args = parse_arguments()
    
    if args.command == 'test':
        print("Running test workflow...")
        return run_script('run_test.sh')
    
    elif args.command == 'batch':
        print("Running batch processing...")
        return run_script('run_batch.sh', genotypes=args.genotypes, animals=args.animals)
    
    elif args.command == 'dry-run':
        print("Running dry run...")
        return run_script('dry_run.sh')
    
    else:
        print("Please specify a command: test, batch, or dry-run")
        return 1

if __name__ == '__main__':
    sys.exit(main())