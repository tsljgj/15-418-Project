#!/usr/bin/env python3
"""
Louvain Algorithm Benchmark Script

This script runs different versions of the Louvain community
detection algorithm and measures performance metrics on heterogeneous systems.
It uses all available resources (P-cores and E-cores) for testing.

Usage:
    python src/checker.py <input_file> [--algorithm ALGORITHM] [--threads THREADS] 
                         [--p-e-ratio P_E_RATIO] [--runs RUNS]

Example:
    python src/checker.py inputs/your_graph_file.txt --algorithm naive,vfc,naive_bl 
           --p-e-ratio 4:12,8:8 --runs 3
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
                      help='Comma-separated list of algorithms to test (sequential,naive,vfc,naive_bl,static,static_bl,vfc_bl) (default: naive,vfc,naive_bl)')
    parser.add_argument('--threads', 
                      help='Comma-separated list of thread counts to test for system-decided cores')
    parser.add_argument('--p-e-ratio', 
                      help='Comma-separated list of P:E core ratios to test (e.g., 4:12,8:8,2:16)')
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
        args.algorithm = 'naive,vfc,naive_bl'
    if args.threads is None and args.p_e_ratio is None:
        # Default P:E ratio of 4:12 (using all 16 cores on an i9-14900K)
        args.p_e_ratio = '4:12'
    
    return args

def run_louvain(executable, input_file, algorithm="sequential", num_threads=None, core_type=None, 
                p_cores=None, e_cores=None):
    """Run the Louvain algorithm and capture output."""
    # Build command
    cmd = [executable, input_file]
    
    # Map algorithm names to command line flags
    algo_flags = {
        "sequential": "-S",
        "naive": "-P",
        "vfc": "-V",
        "naive_bl": "-B",
        "static": "-PS",
        "static_bl": "-PSB",
        "vfc_bl": "-VB"
    }
    
    # Add algorithm flag
    if algorithm in algo_flags:
        cmd.append(algo_flags[algorithm])
    
    # Add core-specific flags
    if algorithm == "sequential":
        if core_type == "p":
            cmd.append('-p')
        elif core_type == "e":
            cmd.append('-e')
    else:
        # For all parallel algorithms
        if p_cores is not None and e_cores is not None:
            cmd.extend(['-pc', str(p_cores), '-ec', str(e_cores)])
        elif num_threads is not None:
            # System decides core allocation
            cmd.extend(['-a', str(num_threads)])
    
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
    
    # Parse configurations
    algorithms = [algo.strip().lower() for algo in args.algorithm.split(',')]
    
    # Parse thread/core configurations
    system_thread_counts = []
    if args.threads:
        system_thread_counts = [int(t) for t in args.threads.split(',')]
        system_thread_counts.sort()
    
    # Parse P:E ratios
    p_e_ratios = []
    if args.p_e_ratio:
        for ratio in args.p_e_ratio.split(','):
            try:
                p, e = map(int, ratio.split(':'))
                p_e_ratios.append((p, e))
            except ValueError:
                print(f"{Fore.RED}Error: Invalid P:E ratio '{ratio}'. Format should be P:E (e.g., 4:12){Style.RESET_ALL}")
                return
    
    # Validate algorithms
    valid_algorithms = ["sequential", "naive", "vfc", "naive_bl", "static", "static_bl", "vfc_bl"]
    algo_display_names = {
        "sequential": "Sequential",
        "naive": "Naive Parallel",
        "vfc": "VFC (Vertex Following + Coloring)",
        "naive_bl": "Naive Parallel + Big.LITTLE",
        "static": "Static Scheduling",
        "static_bl": "Static + Big.LITTLE",
        "vfc_bl": "VFC + Big.LITTLE"
    }
    
    for algo in algorithms:
        if algo not in valid_algorithms:
            print(f"{Fore.RED}Error: Invalid algorithm '{algo}'. Valid options are: {', '.join(valid_algorithms)}{Style.RESET_ALL}")
            return
    
    # Print benchmark configuration
    print_section("CONFIGURATION")
    print(f"Input file:  {Fore.GREEN}{args.input_file}{Style.RESET_ALL}")
    print(f"Algorithms: {Fore.GREEN}{', '.join([algo_display_names.get(algo, algo) for algo in algorithms])}{Style.RESET_ALL}")
    if system_thread_counts:
        print(f"System-decided thread counts: {Fore.GREEN}{', '.join(map(str, system_thread_counts))}{Style.RESET_ALL}")
    if p_e_ratios:
        print(f"P:E core ratios: {Fore.GREEN}{', '.join([f'{p}:{e}' for p, e in p_e_ratios])}{Style.RESET_ALL}")
    print(f"Runs per configuration: {Fore.GREEN}{args.runs}{Style.RESET_ALL}")
    print(f"Executable: {Fore.GREEN}{args.executable}{Style.RESET_ALL}")
    
    # Dictionary to store all results
    all_results = {}
    
    # Run sequential baselines first (always on P and E cores)
    print_section("RUNNING SEQUENTIAL BASELINES")
    sequential_results = []
    for core_type in ['p', 'e']:
        print(f"\nBenchmarking sequential on {Fore.GREEN}{core_type.upper()}-cores{Style.RESET_ALL}:")
        
        runs_runtimes = []
        runs_modularities = []
        
        for i in range(args.runs):
            print(f"  Run {i+1}/{args.runs}...", end="", flush=True)
            stdout, returncode, stderr = run_louvain(args.executable, args.input_file, 
                                                     "sequential", None, core_type)
            
            if returncode == 0:
                runtime, modularity = extract_metrics(stdout)
                runs_runtimes.append(runtime)
                runs_modularities.append(modularity)
                print(f" {Fore.GREEN}Done{Style.RESET_ALL} (Runtime: {runtime:.3f}s, Modularity: {modularity:.6f})")
            else:
                print(f" {Fore.RED}Failed{Style.RESET_ALL}")
        
        # Calculate statistics
        if runs_runtimes:
            avg_runtime = statistics.mean(runs_runtimes)
            avg_modularity = statistics.mean(runs_modularities)
            
            sequential_results.append([
                f"sequential_{core_type}",
                1,
                core_type.upper(),
                "-",
                "-",
                f"{avg_runtime:.3f}",
                f"{avg_modularity:.6f}"
            ])
            
            # Store results
            key = f"sequential_{core_type}"
            all_results[key] = {
                "algorithm": "sequential",
                "threads": 1,
                "core_type": core_type,
                "p_cores": None,
                "e_cores": None,
                "runtime": avg_runtime,
                "modularity": avg_modularity
            }
    
    # Print sequential results
    headers = ["Configuration", "Threads", "Core Type", "P-cores", "E-cores", "Runtime (s)", "Modularity"]
    print_table(headers, sequential_results)
    
    # Run parallel algorithms
    for algorithm in algorithms:
        if algorithm == "sequential":
            continue  # Already done
            
        # Run algorithm but don't print detailed results
        print_section(f"RUNNING {algo_display_names.get(algorithm, algorithm.upper())} ALGORITHM")
        
        # Test with system-decided cores
        if system_thread_counts:
            for threads in system_thread_counts:
                print(f"\nBenchmarking with {Fore.GREEN}{threads}{Style.RESET_ALL} thread(s) (system-decided)...")
                
                runs_runtimes = []
                runs_modularities = []
                
                for i in range(args.runs):
                    stdout, returncode, stderr = run_louvain(args.executable, args.input_file, 
                                                             algorithm, threads)
                    
                    if returncode == 0:
                        runtime, modularity = extract_metrics(stdout)
                        runs_runtimes.append(runtime)
                        runs_modularities.append(modularity)
                    
                # Calculate statistics
                if runs_runtimes:
                    avg_runtime = statistics.mean(runs_runtimes)
                    avg_modularity = statistics.mean(runs_modularities)
                    
                    # Store results
                    key = f"{algorithm}_sys_{threads}"
                    all_results[key] = {
                        "algorithm": algorithm,
                        "threads": threads,
                        "core_type": "system",
                        "p_cores": None,
                        "e_cores": None,
                        "runtime": avg_runtime,
                        "modularity": avg_modularity
                    }
        
        # Test with specific P:E ratios
        if p_e_ratios:
            for p_cores, e_cores in p_e_ratios:
                print(f"\nBenchmarking with {Fore.GREEN}{p_cores} P-cores and {e_cores} E-cores{Style.RESET_ALL}...")
                
                runs_runtimes = []
                runs_modularities = []
                
                for i in range(args.runs):
                    stdout, returncode, stderr = run_louvain(args.executable, args.input_file, 
                                                            algorithm, None, None, 
                                                            p_cores, e_cores)
                    
                    if returncode == 0:
                        runtime, modularity = extract_metrics(stdout)
                        runs_runtimes.append(runtime)
                        runs_modularities.append(modularity)
                
                # Calculate statistics
                if runs_runtimes:
                    avg_runtime = statistics.mean(runs_runtimes)
                    avg_modularity = statistics.mean(runs_modularities)
                    
                    # Store results
                    key = f"{algorithm}_p{p_cores}_e{e_cores}"
                    all_results[key] = {
                        "algorithm": algorithm,
                        "threads": p_cores + e_cores,
                        "core_type": "mixed",
                        "p_cores": p_cores,
                        "e_cores": e_cores,
                        "runtime": avg_runtime,
                        "modularity": avg_modularity
                    }
    
    # Final comparison across all algorithms
    if all_results:
        print_header("PERFORMANCE COMPARISON")
        
        # Prepare comparison table
        comparison_data = []
        
        # Get sequential baselines
        seq_p_runtime = all_results.get("sequential_p", {}).get("runtime", None)
        seq_e_runtime = all_results.get("sequential_e", {}).get("runtime", None)
        
        for key, result in all_results.items():
            if result["algorithm"] == "sequential":
                continue  # Skip sequential in final comparison
                
            # Calculate speedups
            p_speedup = "N/A"
            e_speedup = "N/A"
            
            if seq_p_runtime and result["runtime"] > 0:
                p_speedup = f"{seq_p_runtime / result['runtime']:.2f}x"
            
            if seq_e_runtime and result["runtime"] > 0:
                e_speedup = f"{seq_e_runtime / result['runtime']:.2f}x"
            
            core_config = ""
            if result["core_type"] == "mixed":
                core_config = f"{result['p_cores']}P+{result['e_cores']}E"
            elif result["core_type"] == "system":
                core_config = f"{result['threads']} sys"
            else:
                core_config = result["core_type"].upper()
            
            row = [
                algo_display_names.get(result["algorithm"], result["algorithm"]),
                result["threads"],
                core_config,
                f"{result['runtime']:.3f}",
                f"{result['modularity']:.6f}",
                p_speedup,
                e_speedup
            ]
            
            comparison_data.append(row)
        
        # Sort by algorithm and thread count
        comparison_data.sort(key=lambda x: (x[0], int(x[1])))
        
        # Print comparison table
        headers = ["Algorithm", "Threads", "Core Config", "Runtime (s)", "Modularity", 
                   "Speedup vs P-seq", "Speedup vs E-seq"]
        print_table(headers, comparison_data)

if __name__ == "__main__":
    args = parse_arguments()
    run_benchmarks(args)