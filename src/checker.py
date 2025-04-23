#!/usr/bin/env python3
"""
Louvain Algorithm Benchmark Script

This script runs both sequential and parallel versions of the Louvain community
detection algorithm with varying thread counts and measures performance metrics.

Usage:
    python louvain_benchmark.py <input_file> [--threads THREADS] [--runs RUNS]

Example:
    python louvain_benchmark.py graph.txt --threads 1,2,4,8 --runs 3
"""

import subprocess
import argparse
import re
import os
import time
import statistics
import sys
from colorama import Fore, Style, init

# Initialize colorama
init()

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Benchmark Louvain algorithm implementations.')
    parser.add_argument('input_file', help='Input graph file path')
    parser.add_argument('--threads', 
                      help='Comma-separated list of thread counts to test (default: 1,2,4,8)')
    parser.add_argument('--runs', type=int, default=3, 
                      help='Number of runs for each configuration (default: 3)')
    parser.add_argument('--executable', default='./build/test_louvain',
                      help='Path to the test_louvain executable (default: ./build/test_louvain)')
    
    args = parser.parse_args()
    
    # Set default thread counts if not specified
    if args.threads is None:
        args.threads = '1,2,4,8'
    
    return args

def run_louvain(executable, input_file, parallel=False, num_threads=1):
    """Run the Louvain algorithm and capture output."""
    # Build command
    cmd = [executable, input_file]
    if parallel:
        cmd.extend(['-P', '-n', str(num_threads)])
    else:
        cmd.append('-S')
    
    # Run process and capture output
    try:
        process = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        stdout, stderr = process.communicate(timeout=600)  # 10-minute timeout
        
        # Check for errors
        if process.returncode != 0:
            print(f"{Fore.RED}Error running command: {' '.join(cmd)}{Style.RESET_ALL}")
            print(f"{Fore.RED}Error output:{Style.RESET_ALL}\n{stderr}")
            return None, None, stderr
            
        return stdout, process.returncode, stderr
    except subprocess.TimeoutExpired:
        process.kill()
        print(f"{Fore.RED}Command timed out: {' '.join(cmd)}{Style.RESET_ALL}")
        return None, -1, "Timeout"
    except Exception as e:
        print(f"{Fore.RED}Exception running command: {e}{Style.RESET_ALL}")
        return None, -2, str(e)

def extract_metrics(output):
    """Extract runtime and modularity from command output."""
    if not output:
        return None, None
    
    # Extract total runtime
    runtime_match = re.search(r'Total runtime: (\d+\.\d+) seconds', output)
    runtime = float(runtime_match.group(1)) if runtime_match else None
    
    # Extract final modularity
    modularity_match = re.search(r'Final modularity \(level 0 partition\): ([+-]?\d+\.\d+)', output)
    modularity = float(modularity_match.group(1)) if modularity_match else None
    
    return runtime, modularity

def print_header(header_text):
    """Print a formatted header."""
    width = 80
    print("\n" + "=" * width)
    print(f"{Fore.CYAN}{header_text.center(width)}{Style.RESET_ALL}")
    print("=" * width)

def print_section(section_text):
    """Print a formatted section header."""
    width = 80
    print(f"\n{Fore.YELLOW}{section_text}{Style.RESET_ALL}")
    print("-" * width)

def print_table(headers, data):
    """Print a formatted table without using tabulate."""
    # Find the maximum width needed for each column
    col_widths = [len(h) for h in headers]
    for row in data:
        for i, cell in enumerate(row):
            col_widths[i] = max(col_widths[i], len(str(cell)))
    
    # Print header
    header_row = " | ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
    print("+" + "-" * (len(header_row) + 2) + "+")
    print("| " + header_row + " |")
    print("+" + "-" * (len(header_row) + 2) + "+")
    
    # Print data rows
    for row in data:
        data_row = " | ".join(str(cell).ljust(col_widths[i]) for i, cell in enumerate(row))
        print("| " + data_row + " |")
    
    print("+" + "-" * (len(header_row) + 2) + "+")

def run_benchmarks(args):
    """Run all benchmarks according to provided arguments."""
    print_header("LOUVAIN ALGORITHM BENCHMARK")
    
    # Check if executable exists
    if not os.path.isfile(args.executable):
        print(f"{Fore.RED}Error: Executable not found at {args.executable}{Style.RESET_ALL}")
        print(f"Please build the project first or specify the correct path with --executable.")
        return
    
    # Check if input file exists
    if not os.path.isfile(args.input_file):
        print(f"{Fore.RED}Error: Input file not found at {args.input_file}{Style.RESET_ALL}")
        return
    
    # Parse thread counts
    try:
        thread_counts = [int(t) for t in args.threads.split(',')]
    except ValueError:
        print(f"{Fore.RED}Error: Invalid thread counts. Please provide comma-separated integers.{Style.RESET_ALL}")
        return
    
    # Print benchmark configuration
    print_section("CONFIGURATION")
    print(f"Input file:  {Fore.GREEN}{args.input_file}{Style.RESET_ALL}")
    print(f"Thread counts: {Fore.GREEN}{', '.join(map(str, thread_counts))}{Style.RESET_ALL}")
    print(f"Runs per configuration: {Fore.GREEN}{args.runs}{Style.RESET_ALL}")
    print(f"Executable: {Fore.GREEN}{args.executable}{Style.RESET_ALL}")

    # Run sequential benchmark first
    print_section("RUNNING SEQUENTIAL ALGORITHM")
    seq_runtimes = []
    seq_modularities = []
    
    for i in range(args.runs):
        print(f"Run {i+1}/{args.runs}...", end="", flush=True)
        stdout, returncode, stderr = run_louvain(args.executable, args.input_file, parallel=False)
        
        if returncode == 0:
            runtime, modularity = extract_metrics(stdout)
            seq_runtimes.append(runtime)
            seq_modularities.append(modularity)
            print(f" {Fore.GREEN}Done{Style.RESET_ALL} (Runtime: {runtime:.3f}s, Modularity: {modularity:.6f})")
        else:
            print(f" {Fore.RED}Failed{Style.RESET_ALL}")
    
    # Calculate sequential statistics
    if seq_runtimes:
        avg_seq_runtime = statistics.mean(seq_runtimes)
        avg_seq_modularity = statistics.mean(seq_modularities)
    else:
        print(f"{Fore.RED}Error: All sequential runs failed.{Style.RESET_ALL}")
        return
    
    # Run parallel benchmarks with different thread counts
    print_section("RUNNING PARALLEL ALGORITHM")
    results = []
    
    for threads in thread_counts:
        print(f"\nBenchmarking with {Fore.GREEN}{threads}{Style.RESET_ALL} threads:")
        par_runtimes = []
        par_modularities = []
        
        for i in range(args.runs):
            print(f"  Run {i+1}/{args.runs}...", end="", flush=True)
            stdout, returncode, stderr = run_louvain(args.executable, args.input_file, parallel=True, num_threads=threads)
            
            if returncode == 0:
                runtime, modularity = extract_metrics(stdout)
                par_runtimes.append(runtime)
                par_modularities.append(modularity)
                print(f" {Fore.GREEN}Done{Style.RESET_ALL} (Runtime: {runtime:.3f}s, Modularity: {modularity:.6f})")
            else:
                print(f" {Fore.RED}Failed{Style.RESET_ALL}")
        
        # Calculate statistics for this thread count
        if par_runtimes:
            avg_par_runtime = statistics.mean(par_runtimes)
            avg_par_modularity = statistics.mean(par_modularities)
            speedup = avg_seq_runtime / avg_par_runtime
            
            results.append([
                threads,
                f"{avg_par_runtime:.3f}",
                f"{speedup:.2f}x",
                f"{avg_par_modularity:.6f}",
                f"{(avg_par_modularity - avg_seq_modularity):.6f}"
            ])
        else:
            results.append([threads, "N/A", "N/A", "N/A", "N/A"])
    
    # Print summary table
    print_header("BENCHMARK RESULTS SUMMARY")
    
    # Print sequential baseline
    print(f"Sequential baseline: {Fore.GREEN}{avg_seq_runtime:.3f}s{Style.RESET_ALL}, "
          f"Modularity: {Fore.GREEN}{avg_seq_modularity:.6f}{Style.RESET_ALL}")
    
    # Print parallel results table using our custom table printing function
    headers = ["Threads", "Runtime (s)", "Speedup", "Modularity", "Modularity Diff"]
    print_table(headers, results)
    
    # Print best configuration
    valid_results = [r for r in results if r[1] != "N/A"]
    if valid_results:
        best_speedup_idx = max(range(len(valid_results)), 
                              key=lambda i: float(valid_results[i][2].rstrip('x')))
        best_threads = valid_results[best_speedup_idx][0]
        best_speedup = valid_results[best_speedup_idx][2]
        print(f"\nBest configuration: {Fore.GREEN}{best_threads} threads{Style.RESET_ALL} "
              f"with speedup of {Fore.GREEN}{best_speedup}{Style.RESET_ALL}")
    
    print("\nBenchmark completed successfully.")

if __name__ == "__main__":
    args = parse_arguments()
    run_benchmarks(args)