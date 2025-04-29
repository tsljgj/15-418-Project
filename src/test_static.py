import subprocess
import time
import os
import csv
import sys
import re
from itertools import product

def run_test(command):
    """Run a test and capture its execution time and modularity"""
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    
    # Parse the output to get modularity and runtime
    runtime_match = re.search(r"Total runtime: (\d+\.\d+) seconds", result.stdout)
    modularity_match = re.search(r"Final modularity \(level 0 partition\): (\d+\.\d+)", result.stdout)
    
    reported_runtime = float(runtime_match.group(1)) if runtime_match else None
    modularity = float(modularity_match.group(1)) if modularity_match else None
    
    return {
        "reported_runtime": reported_runtime,
        "modularity": modularity
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python test_static.py <graph_file>")
        sys.exit(1)
    
    graph_file = sys.argv[1]
    base_path = "build\\test_louvain.exe"
    graph_path = graph_file
    file_name = os.path.basename(graph_file)
    
    # Create results directory if it doesn't exist
    os.makedirs("results", exist_ok=True)
    
    # Prepare CSV file for results
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    results_file = f"results\\speed_ratio_test_{file_name}_{timestamp}.csv"
    
    with open(results_file, 'w', newline='') as csvfile:
        fieldnames = ['test_type', 'p_cores', 'e_cores', 'speed_ratio', 
                      'speedup', 'modularity']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        # Run baseline sequential tests
        print("=== Running P-core baseline ===")
        p_baseline = run_test(f"{base_path} {graph_path} -S -p")
        
        print("\n=== Running E-core baseline ===")
        e_baseline = run_test(f"{base_path} {graph_path} -S -e")
        
        # Store E-core baseline runtime for speedup calculations
        e_baseline_runtime = e_baseline['reported_runtime']
        
        # Calculate and save P-core speedup
        p_speedup = e_baseline_runtime / p_baseline['reported_runtime'] if p_baseline['reported_runtime'] else None
        
        writer.writerow({
            'test_type': 'sequential_p',
            'p_cores': 1, 'e_cores': 0, 'speed_ratio': 'N/A',
            'speedup': p_speedup,
            'modularity': p_baseline['modularity']
        })
        
        writer.writerow({
            'test_type': 'sequential_e',
            'p_cores': 0, 'e_cores': 1, 'speed_ratio': 'N/A',
            'speedup': 1.000,  # By definition
            'modularity': e_baseline['modularity']
        })
        
        # Calculate actual observed speed ratio
        if e_baseline['reported_runtime'] and p_baseline['reported_runtime']:
            actual_ratio = e_baseline['reported_runtime'] / p_baseline['reported_runtime']
            print(f"\nObserved P:E speed ratio = {actual_ratio:.3f}")
            
            writer.writerow({
                'test_type': 'observed_ratio',
                'p_cores': 'N/A', 'e_cores': 'N/A', 
                'speed_ratio': actual_ratio,
                'speedup': 'N/A',
                'modularity': 'N/A'
            })
            
            csvfile.flush()
        
        # Run parallel tests with various configurations
        p_cores_range = [1, 2, 4]  # As requested
        e_cores_range = [1, 2, 4, 8, 16]  # As requested
        speed_ratios = [1.5, 2.0, 2.5, 3.0, 3.5]  # As requested
        
        # Test all combinations of core counts with different speed ratios
        combinations = list(product(p_cores_range, e_cores_range, speed_ratios))
        total_tests = len(combinations)
        
        for i, (p_count, e_count, ratio) in enumerate(combinations):
            print(f"Test {i+1}/{total_tests}: P={p_count}, E={e_count}, R={ratio}")
            result = run_test(f"{base_path} {graph_path} -PSB -pc {p_count} -ec {e_count} -sr {ratio}")
            
            # Calculate speedup compared to E-core baseline
            speedup = e_baseline_runtime / result['reported_runtime'] if result['reported_runtime'] else None
            
            # Format speedup to preserve 3 decimal places
            if speedup:
                speedup_str = f"{speedup:.3f}"
            else:
                speedup_str = None
                
            print(f"Speedup: {speedup_str}x, Modularity: {result['modularity']}")
            
            writer.writerow({
                'test_type': 'parallel_static_bl',
                'p_cores': p_count, 
                'e_cores': e_count, 
                'speed_ratio': ratio,
                'speedup': speedup_str,
                'modularity': result['modularity']
            })
            csvfile.flush()
    
    print(f"\nAll tests completed. Results saved to {results_file}")
    
    # Find best performance
    try:
        import pandas as pd
        df = pd.read_csv(results_file)
        parallel_results = df[df['test_type'] == 'parallel_static_bl'].copy()
        
        if not parallel_results.empty:
            # Convert speedup from string to float for proper sorting
            parallel_results['speedup'] = parallel_results['speedup'].astype(float)
            
            best_config = parallel_results.loc[parallel_results['speedup'].idxmax()]
            print(f"\nBest configuration:")
            print(f"  P-cores: {best_config['p_cores']}, E-cores: {best_config['e_cores']}")
            print(f"  Speed ratio: {best_config['speed_ratio']}")
            print(f"  Speedup: {best_config['speedup']}x")
            print(f"  Modularity: {best_config['modularity']}")
            
            # Print top 5 configurations
            top5 = parallel_results.nlargest(5, 'speedup')
            print("\nTop 5 configurations:")
            for idx, row in top5.iterrows():
                print(f"  P-cores: {row['p_cores']}, E-cores: {row['e_cores']}, "
                      f"Ratio: {row['speed_ratio']}, Speedup: {row['speedup']}x")
    except ImportError:
        print("Pandas not installed. Cannot analyze results automatically.")
        print("Please analyze the CSV file manually to find the best configuration.")

if __name__ == "__main__":
    main()