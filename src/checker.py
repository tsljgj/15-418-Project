#!/usr/bin/env python3
"""
Louvain Algorithm Benchmark Script

This script runs different versions of the Louvain community
detection algorithm with varying thread counts and measures performance metrics.

Usage:
    python src/checker.py <input_file> [--algorithm ALGORITHM] [--threads THREADS] [--runs RUNS]
    [--core-type CORE_TYPE]

Example:
    python src/checker.py inputs/community_graph_5e5.txt --algorithm sequential,naive,vfc --threads 1,2,4,8 --runs 1 --core-type p,e
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
    parser.add_argument('--core-type', 
                      help='Comma-separated list of core types to test for sequential algorithm (p,e) (default: p)')

    # Automatically detect if .exe exists
    if os.path.isfile('./build/Release/test_louvain.exe'):
        exe_path = './build/Release/test_louvain.exe'
    elif os.path.isfile('./build/test_louvain.exe'):
        exe_path = './build/test_louvain.exe'
    elif os.path.isfile('./build/test_louvain'):
        exe_path = './build/test_louvain'
    else:
        exe_path = './build/Release/test_louvain'
        
    parser.add_argument('--executable', default=exe_path,
                      help='Path to the test_louvain executable (auto-detects location)')
    
    args = parser.parse_args()

    # Set defaults if not provided
    if args.algorithm is None:
        args.algorithm = 'sequential,naive,vfc'
    if args.threads is None:
        args.threads = '1,2,4,8'
    if args.core_type is None:
        args.core_type = 'p'
    
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
    
    # Parse algorithms
    algorithms = [algo.strip().lower() for algo in args.algorithm.split(',')]
    valid_algorithms = ["sequential", "naive", "vfc"]
    for algo in algorithms:
        if algo not in valid_algorithms:
            print(f"{Fore.RED}Error: Invalid algorithm '{algo}'. Valid options are: {', '.join(valid_algorithms)}{Style.RESET_ALL}")
            return
    
    # Parse core types
    core_types = [ct.strip().lower() for ct in args.core_type.split(',')]
    valid_core_types = ["p", "e"]
    for ct in core_types:
        if ct not in valid_core_types:
            print(f"{Fore.RED}Error: Invalid core type '{ct}'. Valid options are: {', '.join(valid_core_types)}{Style.RESET_ALL}")
            return
    
    # Print benchmark configuration
    print_section("CONFIGURATION")
    print(f"Input file:  {Fore.GREEN}{args.input_file}{Style.RESET_ALL}")
    print(f"Algorithms: {Fore.GREEN}{', '.join(algorithms)}{Style.RESET_ALL}")
    print(f"Thread counts: {Fore.GREEN}{', '.join(map(str, thread_counts))}{Style.RESET_ALL}")
    print(f"Core types (for sequential): {Fore.GREEN}{', '.join(core_types)}{Style.RESET_ALL}")
    print(f"Runs per configuration: {Fore.GREEN}{args.runs}{Style.RESET_ALL}")
    print(f"Executable: {Fore.GREEN}{args.executable}{Style.RESET_ALL}")

    # Dictionary to store all results for final comparison
    all_results = {}
    
    # Run benchmarks for each algorithm
    for algorithm in algorithms:
        print_section(f"RUNNING {algorithm.upper()} ALGORITHM")
        
        # For sequential, we only need to run it once with each core type
        if algorithm == "sequential":
            thread_counts_to_use = [1]
            core_types_to_use = core_types
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
    
    # Final comparison across all algorithms
    if all_results:
        print_header("COMPARISON ACROSS ALL ALGORITHMS AND CORE TYPES")
        
        # Prepare comparison table
        comparison_data = []
        for key, result in all_results.items():
            algorithm = result["algorithm"]
            threads = result["threads"]
            core_type = result["core_type"] if result["core_type"] else "N/A"
            runtime = result["runtime"]
            modularity = result["modularity"]
            
            comparison_data.append([
                algorithm,
                threads,
                core_type.upper() if core_type != "N/A" else core_type,
                f"{runtime:.3f}",
                f"{modularity:.6f}"
            ])
        
        # Sort comparison_data by algorithm name, thread count, and core type
        comparison_data.sort(key=lambda x: (x[0], x[1], x[2]))
        
        # Print comparison table
        headers = ["Algorithm", "Threads", "Core Type", "Runtime (s)", "Modularity"]
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
        
        # If sequential is part of the benchmark with both P and E cores, make a direct comparison
        if "sequential" in algorithms and len(core_types) > 1:
            print("\nDirect P-core vs E-core comparison (Sequential algorithm):")
            
            p_core_key = next((k for k, v in all_results.items() 
                            if v["algorithm"] == "sequential" and v["core_type"] == "p"), None)
            e_core_key = next((k for k, v in all_results.items() 
                            if v["algorithm"] == "sequential" and v["core_type"] == "e"), None)
            
            if p_core_key and e_core_key:
                p_core_runtime = all_results[p_core_key]["runtime"]
                e_core_runtime = all_results[e_core_key]["runtime"]
                p_core_modularity = all_results[p_core_key]["modularity"]
                e_core_modularity = all_results[e_core_key]["modularity"]
                
                runtime_diff = abs(p_core_runtime - e_core_runtime)
                runtime_ratio = e_core_runtime / p_core_runtime if p_core_runtime > 0