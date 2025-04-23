#!/usr/bin/env python3
"""
Louvain Algorithm Benchmark Script

This script runs different versions of the Louvain community
detection algorithm with varying thread counts and measures performance metrics.

Usage:
    python src/checker.py <input_file> [--algorithm ALGORITHM] [--threads THREADS] [--runs RUNS]

Example:
    python src/checker.py inputs/community_graph_5e5.txt --algorithm sequential,naive,vfc --threads 1,2,4,8 --runs 1
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
    parser.add_argument('--algorithm', 
                      help='Comma-separated list of algorithms to test (sequential,naive,vfc) (default: sequential,naive,vfc)')
    parser.add_argument('--threads', 
                      help='Comma-separated list of thread counts to test (default: 1,2,4,8)')
    parser.add_argument('--runs', type=int, default=1, 
                      help='Number of runs for each configuration (default: 1)')
    parser.add_argument('--executable', default='./build/test_louvain',
                      help='Path to the test_louvain executable (default: ./build/test_louvain)')
    
    args = parser.parse_args()
    
    # Set default algorithms if not specified
    if args.algorithm is None:
        args.algorithm = 'sequential,naive,vfc'
    
    # Set default thread counts if not specified
    if args.threads is None:
        args.threads = '1,2,4,8'
    
    return args

def run_louvain(executable, input_file, algorithm="sequential", num_threads=1):
    """Run the Louvain algorithm and capture output."""
    # Build command
    cmd = [executable, input_file]
    if algorithm == "sequential":
        cmd.append('-S')
    elif algorithm == "naive":
        cmd.extend(['-P', '-n', str(num_threads)])
    elif algorithm == "vfc":
        cmd.extend(['-V', '-n', str(num_threads)])
    
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
    
    # Parse algorithms
    algorithms = [algo.strip().lower() for algo in args.algorithm.split(',')]
    valid_algorithms = ["sequential", "naive", "vfc"]
    for algo in algorithms:
        if algo not in valid_algorithms:
            print(f"{Fore.RED}Error: Invalid algorithm '{algo}'. Valid options are: {', '.join(valid_algorithms)}{Style.RESET_ALL}")
            return
    
    # Print benchmark configuration
    print_section("CONFIGURATION")
    print(f"Input file:  {Fore.GREEN}{args.input_file}{Style.RESET_ALL}")
    print(f"Algorithms: {Fore.GREEN}{', '.join(algorithms)}{Style.RESET_ALL}")
    print(f"Thread counts: {Fore.GREEN}{', '.join(map(str, thread_counts))}{Style.RESET_ALL}")
    print(f"Runs per configuration: {Fore.GREEN}{args.runs}{Style.RESET_ALL}")
    print(f"Executable: {Fore.GREEN}{args.executable}{Style.RESET_ALL}")

    # Dictionary to store all results for final comparison
    all_results = {}
    
    # Run benchmarks for each algorithm
    for algorithm in algorithms:
        print_section(f"RUNNING {algorithm.upper()} ALGORITHM")
        
        # For sequential, we only need to run it once
        if algorithm == "sequential":
            thread_counts_to_use = [1]
        else:
            thread_counts_to_use = thread_counts
        
        algorithm_results = []
        
        for threads in thread_counts_to_use:
            print(f"\nBenchmarking with {Fore.GREEN}{threads}{Style.RESET_ALL} thread(s):")
            runs_runtimes = []
            runs_modularities = []
            
            for i in range(args.runs):
                print(f"  Run {i+1}/{args.runs}...", end="", flush=True)
                stdout, returncode, stderr = run_louvain(args.executable, args.input_file, algorithm, threads)
                
                if returncode == 0:
                    runtime, modularity = extract_metrics(stdout)
                    runs_runtimes.append(runtime)
                    runs_modularities.append(modularity)
                    print(f" {Fore.GREEN}Done{Style.RESET_ALL} (Runtime: {runtime:.3f}s, Modularity: {modularity:.6f})")
                else:
                    print(f" {Fore.RED}Failed{Style.RESET_ALL}")
            
            # Calculate statistics for this thread count
            if runs_runtimes:
                avg_runtime = statistics.mean(runs_runtimes)
                avg_modularity = statistics.mean(runs_modularities)
                
                algorithm_results.append([
                    threads,
                    f"{avg_runtime:.3f}",
                    f"{avg_modularity:.6f}"
                ])
                
                # Store for final comparison
                key = f"{algorithm}_{threads}"
                all_results[key] = {
                    "algorithm": algorithm,
                    "threads": threads,
                    "runtime": avg_runtime,
                    "modularity": avg_modularity
                }
            else:
                algorithm_results.append([threads, "N/A", "N/A"])
        
        # Print results table for this algorithm
        headers = ["Threads", "Runtime (s)", "Modularity"]
        print_table(headers, algorithm_results)
    
    # Final comparison across all algorithms
    if all_results:
        print_header("COMPARISON ACROSS ALL ALGORITHMS")
        
        # Find sequential runtime and modularity for baseline
        seq_runtime = None
        seq_modularity = None
        for key, result in all_results.items():
            if result["algorithm"] == "sequential":
                seq_runtime = result["runtime"]
                seq_modularity = result["modularity"]
                break
        
        # Print sequential baseline info
        if seq_runtime is not None:
            print(f"Sequential baseline: Runtime = {Fore.GREEN}{seq_runtime:.3f}s{Style.RESET_ALL}, "
                  f"Modularity = {Fore.GREEN}{seq_modularity:.6f}{Style.RESET_ALL}")
            print()
        
        # Find the best modularity and runtime across all parallel algorithms
        parallel_results = {k: v for k, v in all_results.items() if v["algorithm"] != "sequential"}
        
        if not parallel_results:
            print("No parallel algorithm results to compare.")
            return
            
        best_modularity = max(result["modularity"] for result in parallel_results.values())
        best_runtime = min(result["runtime"] for result in parallel_results.values())
        
        # Prepare comparison table (excluding sequential)
        comparison_data = []
        for key, result in sorted(parallel_results.items()):
            algorithm = result["algorithm"]
            threads = result["threads"]
            runtime = result["runtime"]
            modularity = result["modularity"]
            
            # Calculate speedup relative to sequential
            if seq_runtime is not None:
                speedup = f"{seq_runtime / runtime:.2f}x"
            else:
                speedup = "N/A"
            
            # Calculate modularity relative to best
            modularity_rel = f"{modularity / best_modularity:.4f}"
            
            comparison_data.append([
                algorithm,
                threads,
                f"{runtime:.3f}",
                speedup,
                f"{modularity:.6f}",
                modularity_rel
            ])
        
        # Print comparison table
        headers = ["Algorithm", "Threads", "Runtime (s)", "Speedup", "Modularity", "Modularity (rel)"]
        print_table(headers, comparison_data)
        
        # Find the best configurations
        best_modularity_config = max(all_results.items(), key=lambda x: x[1]["modularity"])
        best_runtime_config = min(all_results.items(), key=lambda x: x[1]["runtime"])
        
        print(f"\nBest modularity: {Fore.GREEN}{best_modularity_config[1]['modularity']:.6f}{Style.RESET_ALL} "
              f"with {Fore.GREEN}{best_modularity_config[1]['algorithm']}{Style.RESET_ALL} "
              f"using {Fore.GREEN}{best_modularity_config[1]['threads']}{Style.RESET_ALL} threads")
        
        print(f"Best runtime: {Fore.GREEN}{best_runtime_config[1]['runtime']:.3f}s{Style.RESET_ALL} "
              f"with {Fore.GREEN}{best_runtime_config[1]['algorithm']}{Style.RESET_ALL} "
              f"using {Fore.GREEN}{best_runtime_config[1]['threads']}{Style.RESET_ALL} threads")
    
    print("\nBenchmark completed successfully.")

if __name__ == "__main__":
    args = parse_arguments()
    run_benchmarks(args)