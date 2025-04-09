#include "louvain_seq.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <unordered_set>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>\n";
        return EXIT_FAILURE;
    }
    std::string input_filename = argv[1];
    
    // Read the graph.
    Graph g = readGraph(input_filename);
    std::cout << "Graph loaded: " << g.n << " nodes, total weighted m = " << g.m << "\n";
    
    // Run hierarchical Louvain algorithm and time its execution.
    auto start_time = std::chrono::high_resolution_clock::now();
    
    Hierarchy H;
    louvainHierarchical(g, H);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    
    std::cout << "\n--- Hierarchical Decomposition Summary ---\n";
    for (size_t level = 0; level < H.partitions.size(); ++level) {
        const std::vector<int> &part = H.partitions[level];
        std::unordered_set<int> comms(part.begin(), part.end());
        std::cout << "Level " << level << " -> " << comms.size() << " communities\n";
    }
    
    std::cout << "\nTotal runtime: " << elapsed_seconds.count() << " seconds\n";
    
    // Also output final modularity computed on the original graph.
    double Q = computeModularity(g, H.partitions[0]);
    std::cout << "Final modularity (level 0 partition): " << Q << "\n";
    
    return EXIT_SUCCESS;
}