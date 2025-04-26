#!/usr/bin/env python3
"""
Louvain Algorithm Benchmark Script

This script runs different versions of the Louvain community
detection algorithm with varying thread counts and measures performance metrics.
It automatically establishes sequential baselines on both P-cores and E-cores
for comparison.

Usage:
    python src/checker.py <input_file> [--algorithm ALGORITHM] [--threads THREADS] [--runs RUNS]

Example:
    python src/checker.py inputs/community_graph_5e5.txt --algorithm naive,vfc --threads 1,2,4,8 --runs 1
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
                      help='Comma-separated list of algorithms to test (naive,vfc) (default: naive,vfc)')
    parser.add_argument('--threads', 
                      help='Comma-separated list of thread counts to test (default: 1,2,4,8)')
    parser.add_argument('--runs', type=int, default=1, 
                      help='Number of runs for each configuration (default: 1)')

    # Automatically detect if .exe exists
    if os.path.isfile('./build/Release/test_louvain.exe'):
        exe_path = './build/Release/test_louvain.exe'
    elif os.path.isfile('./build/test_louvain.exe'):
        exe_path = './build/test_louvain.exe'
    elif os.path.isfile('./build/test_louvain'):
        exe_path = './build/test_louvain'
    else:
        exe_path = './build/test_louvain.exe'  # Default to build directory
        
    parser.add_argument('--executable', default=exe_path,
                      help='Path to the test_louvain executable (auto-detects location)')
    
    args = parser.parse_args()

    # Set defaults if not provided
    if args.algorithm is None:
        args.algorithm = 'naive,vfc'
    if args.threads is None:
        args.threads = '1,2,4,8'
    
    return args

def run_louvain(executable, input_file, algorithm="sequential", num_threads=1, core_type=None):
    """Run the Louvain algorithm and capture output."""
    # Build command
    cmd = [executable, input_file]
    if algorithm == "sequential":
        cmd.append('-S')
        if core_type == "p":
            cmd.append('-p')
        elif core_type == "e":
            cmd.append('-e')
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
        # Sort thread counts to ensure they appear in ascending order
        thread_counts.sort()
    except ValueError:
        print(f"{Fore.RED}Error: Invalid thread counts. Please provide comma-separated integers.{Style.RESET_ALL}")
        return
    
    # Parse algorithms - always include 'sequential' internally
    requested_algorithms = [algo.strip().lower() for algo in args.algorithm.split(',')]
    algorithms = ['sequential']  # Always run sequential first for baselines
    algorithms.extend([algo for algo in requested_algorithms if algo != 'sequential'])
    
    valid_algorithms = ["sequential", "naive", "vfc"]
    for algo in algorithms:
        if algo not in valid_algorithms:
            print(f"{Fore.RED}Error: Invalid algorithm '{algo}'. Valid options are: {', '.join(valid_algorithms)}{Style.RESET_ALL}")
            return
    
    # Core types - always run both P and E cores for sequential
    core_types = ["p", "e"]
    
    # Print benchmark configuration
    print_section("CONFIGURATION")
    print(f"Input file:  {Fore.GREEN}{args.input_file}{Style.RESET_ALL}")
    print(f"Requested algorithms: {Fore.GREEN}{', '.join(requested_algorithms)}{Style.RESET_ALL}")
    print(f"Thread counts: {Fore.GREEN}{', '.join(map(str, thread_counts))}{Style.RESET_ALL}")
    print(f"Sequential baselines: {Fore.GREEN}P-core and E-core{Style.RESET_ALL}")
    print(f"Runs per configuration: {Fore.GREEN}{args.runs}{Style.RESET_ALL}")
    print(f"Executable: {Fore.GREEN}{args.executable}{Style.RESET_ALL}")

    # Dictionary to store all results for final comparison
    all_results = {}
    
    # Variables to store sequential baseline results
    sequential_p_runtime = None
    sequential_e_runtime = None
    sequential_p_modularity = None
    sequential_e_modularity = None
    
    # Run benchmarks for each algorithm
    for algorithm in algorithms:
        print_section(f"RUNNING {algorithm.upper()} ALGORITHM")
        
        # For sequential, we need to run it with both P and E cores
        if algorithm == "sequential":
            thread_counts_to_use = [1]  # Sequential always uses 1 thread
            core_types_to_use = core_types  # Test both P and E cores
        else:
            thread_counts_to_use = thread_counts
            core_types_to_use = [None]  # No core type for parallel algorithms
        
        algorithm_results = []
        
        for core_type in core_types_to_use:
            for threads in thread_counts_to_use:
                print_info = f"\nBenchmarking with {Fore.GREEN}{threads}{Style.RESET_ALL} thread(s)"
                if core_type:
                    print_info += f" on {Fore.GREEN}{core_type.upper()}-cores{Style.RESET_ALL}"
                print(print_info + ":")
                
                runs_runtimes = []
                runs_modularities = []
                
                for i in range(args.runs):
                    print(f"  Run {i+1}/{args.runs}...", end="", flush=True)
                    stdout, returncode, stderr = run_louvain(args.executable, args.input_file, algorithm, threads, core_type)
                    
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
                    
                    # Add core type information for sequential algorithm
                    if core_type:
                        result_row = [
                            threads,
                            core_type.upper(),
                            f"{avg_runtime:.3f}",
                            f"{avg_modularity:.6f}"
                        ]
                    else:
                        result_row = [
                            threads,
                            "N/A",
                            f"{avg_runtime:.3f}",
                            f"{avg_modularity:.6f}"
                        ]
                    
                    algorithm_results.append(result_row)
                    
                    # Store for final comparison
                    key = f"{algorithm}_{threads}_{core_type}" if core_type else f"{algorithm}_{threads}"
                    all_results[key] = {
                        "algorithm": algorithm,
                        "threads": threads,
                        "core_type": core_type,
                        "runtime": avg_runtime,
                        "modularity": avg_modularity
                    }
                    
                    # Store sequential baselines
                    if algorithm == "sequential" and threads == 1:
                        if core_type == "p":
                            sequential_p_runtime = avg_runtime
                            sequential_p_modularity = avg_modularity
                        elif core_type == "e":
                            sequential_e_runtime = avg_runtime
                            sequential_e_modularity = avg_modularity
                else:
                    if core_type:
                        algorithm_results.append([threads, core_type.upper(), "N/A", "N/A"])
                    else:
                        algorithm_results.append([threads, "N/A", "N/A", "N/A"])
        
        # Sort algorithm_results by thread count before printing
        algorithm_results.sort(key=lambda x: (x[0], x[1]))
        
        # Print results table for this algorithm
        if algorithm == "sequential":
            headers = ["Threads", "Core Type", "Runtime (s)", "Modularity"]
        else:
            headers = ["Threads", "Core Type", "Runtime (s)", "Modularity"]
        print_table(headers, algorithm_results)
        
        # For sequential, print a direct comparison between P and E cores
        if algorithm == "sequential" and sequential_p_runtime and sequential_e_runtime:
            print("\nDirect P-core vs E-core comparison (Sequential algorithm):")
            
            runtime_diff = abs(sequential_p_runtime - sequential_e_runtime)
            runtime_ratio = sequential_p_runtime / sequential_e_runtime if sequential_e_runtime > 0 else float('inf')
            
            print(f"P-core runtime: {Fore.GREEN}{sequential_p_runtime:.3f}s{Style.RESET_ALL}")
            print(f"E-core runtime: {Fore.GREEN}{sequential_e_runtime:.3f}s{Style.RESET_ALL}")
            print(f"Runtime difference: {Fore.GREEN}{runtime_diff:.3f}s{Style.RESET_ALL}")
            print(f"P-core to E-core ratio: {Fore.GREEN}{runtime_ratio:.3f}x{Style.RESET_ALL}")
            
            if sequential_p_modularity != sequential_e_modularity:
                print(f"Note: Modularity values differ: P-core = {sequential_p_modularity:.6f}, "
                      f"E-core = {sequential_e_modularity:.6f}")
    
    # Skip output for sequential if it wasn't requested, to focus on the algorithms that were requested
    if 'sequential' not in requested_algorithms:
        # Remove sequential results from all_results for final comparison table
        all_results = {k: v for k, v in all_results.items() if v['algorithm'] != 'sequential'}
    
    # Final comparison across all algorithms
    if all_results:
        print_header("PERFORMANCE COMPARISON")
        
        # Prepare comparison table
        comparison_data = []
        for key, result in all_results.items():
            algorithm = result["algorithm"]
            # Skip sequential baselines if they weren't explicitly requested
            if algorithm == "sequential" and algorithm not in requested_algorithms:
                continue
                
            threads = result["threads"]
            core_type = result["core_type"] if result["core_type"] else "N/A"
            runtime = result["runtime"]
            modularity = result["modularity"]
            
            # Calculate speedups if we have sequential baselines
            p_speedup = "N/A"
            e_speedup = "N/A"
            
            if sequential_p_runtime and runtime > 0:
                p_speedup = f"{sequential_p_runtime / runtime:.2f}x"
                
            if sequential_e_runtime and runtime > 0:
                e_speedup = f"{sequential_e_runtime / runtime:.2f}x"
                
            row = [
                algorithm,
                threads,
                core_type.upper() if core_type != "N/A" else core_type,
                f"{runtime:.3f}",
                f"{modularity:.6f}",
                p_speedup,
                e_speedup
            ]
            
            comparison_data.append(row)
        
        # Sort by algorithm name, thread count, and core type
        comparison_data.sort(key=lambda x: (x[0], int(x[1]), x[2]))
        
        # Print comparison table
        headers = ["Algorithm", "Threads", "Core Type", "Runtime (s)", "Modularity", 
                   "Speedup vs P-core", "Speedup vs E-core"]
        print_table(headers, comparison_data)
        
        # Find the best configurations
        best_modularity_config = max(all_results.items(), key=lambda x: x[1]["modularity"])
        best_runtime_config = min(all_results.items(), key=lambda x: x[1]["runtime"])
        
        print(f"\nBest modularity: {Fore.GREEN}{best_modularity_config[1]['modularity']:.6f}{Style.RESET_ALL} "
              f"with {Fore.GREEN}{best_modularity_config[1]['algorithm']}{Style.RESET_ALL} "
              f"using {Fore.GREEN}{best_modularity_config[1]['threads']}{Style.RESET_ALL} threads")
        
        if best_modularity_config[1]['core_type']:
            print(f"  on {Fore.GREEN}{best_modularity_config[1]['core_type'].upper()}-cores{Style.RESET_ALL}")
        
        print(f"Best runtime: {Fore.GREEN}{best_runtime_config[1]['runtime']:.3f}s{Style.RESET_ALL} "
              f"with {Fore.GREEN}{best_runtime_config[1]['algorithm']}{Style.RESET_ALL} "
              f"using {Fore.GREEN}{best_runtime_config[1]['threads']}{Style.RESET_ALL} threads")
        
        if best_runtime_config[1]['core_type']:
            print(f"  on {Fore.GREEN}{best_runtime_config[1]['core_type'].upper()}-cores{Style.RESET_ALL}")
        
        # Speedup summary
        print_section("SPEEDUP SUMMARY")
        
        max_p_speedup = 0
        max_e_speedup = 0
        max_p_config = None
        max_e_config = None
        
        for key, result in all_results.items():
            # Skip sequential results
            if result["algorithm"] == "sequential" and result["threads"] == 1:
                continue
                
            runtime = result["runtime"]
            if runtime > 0:
                if sequential_p_runtime:
                    p_speedup = sequential_p_runtime / runtime
                    if p_speedup > max_p_speedup:
                        max_p_speedup = p_speedup
                        max_p_config = result
                        
                if sequential_e_runtime:
                    e_speedup = sequential_e_runtime / runtime
                    if e_speedup > max_e_speedup:
                        max_e_speedup = e_speedup
                        max_e_config = result
        
        if max_p_config:
            print(f"Maximum speedup vs P-core sequential: {Fore.GREEN}{max_p_speedup:.2f}x{Style.RESET_ALL}")
            print(f"  with {Fore.GREEN}{max_p_config['algorithm']}{Style.RESET_ALL} "
                  f"using {Fore.GREEN}{max_p_config['threads']}{Style.RESET_ALL} threads")
            if max_p_config['core_type']:
                print(f"  on {Fore.GREEN}{max_p_config['core_type'].upper()}-cores{Style.RESET_ALL}")
                
        if max_e_config:
            print(f"Maximum speedup vs E-core sequential: {Fore.GREEN}{max_e_speedup:.2f}x{Style.RESET_ALL}")
            print(f"  with {Fore.GREEN}{max_e_config['algorithm']}{Style.RESET_ALL} "
                  f"using {Fore.GREEN}{max_e_config['threads']}{Style.RESET_ALL} threads")
            if max_e_config['core_type']:
                print(f"  on {Fore.GREEN}{max_e_config['core_type'].upper()}-cores{Style.RESET_ALL}")

if __name__ == "__main__":
    args = parse_arguments()
    run_benchmarks(args)