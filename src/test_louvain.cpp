#include "louvain_seq.h"
#include "louvain_parallel.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <unordered_set>
#include <cstring>

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " <input_file> [-S|-P [-n num_threads]]\n";
    std::cerr << "Options:\n";
    std::cerr << "  -S               Run sequential Louvain algorithm (default)\n";
    std::cerr << "  -P               Run parallel Louvain algorithm with graph partitioning\n";
    std::cerr << "  -n num_threads   Number of threads/partitions to use (default: 1)\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    std::string input_filename = argv[1];
    bool useParallel = false;
    int numThreads = 1;
    
    // Parse command line arguments
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-S") == 0) {
            useParallel = false;
        }
        else if (strcmp(argv[i], "-P") == 0) {
            useParallel = true;
        }
        else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            numThreads = std::atoi(argv[++i]);
            if (numThreads <= 0) {
                std::cerr << "Error: Number of threads must be positive\n";
                return EXIT_FAILURE;
            }
        }
        else {
            printUsage(argv[0]);
            return EXIT_FAILURE;
        }
    }
    
    // Read the graph
    Graph g = readGraph(input_filename);
    std::cout << "Graph loaded: " << g.n << " nodes, total weighted m = " << g.m << "\n";
    
    // Run the appropriate algorithm and time its execution
    auto start_time = std::chrono::high_resolution_clock::now();
    
    Hierarchy H;
    if (useParallel) {
        std::cout << "Running parallel Louvain algorithm with " << numThreads << " threads\n";
        louvainParallel(g, H, numThreads);
    } else {
        std::cout << "Running sequential Louvain algorithm\n";
        louvainHierarchical(g, H);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    
    std::cout << "\n--- Hierarchical Decomposition Summary ---\n";
    for (size_t level = 0; level < H.partitions.size(); ++level) {
        const std::vector<int> &part = H.partitions[level];
        std::unordered_set<int> comms(part.begin(), part.end());
        std::cout << "Level " << level << " -> " << comms.size() << " communities\n";
    }
    
    std::cout << "\nTotal runtime: " << elapsed_seconds.count() << " seconds\n";
    
    // Also output final modularity computed on the original graph
    double Q = computeModularity(g, H.partitions[0]);
    std::cout << "Final modularity (level 0 partition): " << Q << "\n";
    
    return EXIT_SUCCESS;
}