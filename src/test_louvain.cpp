#include "louvain_seq.h"
#include "louvain_parallel.h"
#include "louvain_parallel_vfc.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <unordered_set>
#include <cstring>

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " <input_file> [-S|-P|-V] [options]\n";
    std::cerr << "Options:\n";
    std::cerr << "  -S               Run sequential Louvain algorithm (default)\n";
    std::cerr << "  -P               Run naive parallel Louvain algorithm\n";
    std::cerr << "  -V               Run parallel Louvain algorithm with Vertex Following and Coloring\n";
    std::cerr << "  For sequential algorithm:\n";
    std::cerr << "    -p             Run on P-cores (performance cores)\n";
    std::cerr << "    -e             Run on E-cores (efficiency cores)\n";
    std::cerr << "    -a             Use any available cores (system decides, default)\n";
    std::cerr << "  For parallel algorithms:\n";
    std::cerr << "    -pc num        Number of P-cores to use\n";
    std::cerr << "    -ec num        Number of E-cores to use\n";
    std::cerr << "    -a num         Use 'num' threads on any available cores (system decides)\n";
    std::cerr << "    -n num         DEPRECATED: Same as -a num\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    std::string input_filename = argv[1];
    enum AlgorithmType { SEQUENTIAL, NAIVE_PARALLEL, PARALLEL_VFC };
    AlgorithmType algorithm = SEQUENTIAL;
    int numThreads = 1;
    
    // New variables for P-core and E-core counts
    int pCoreCount = 0;
    int eCoreCount = 0;
    bool useSystemCores = true; // Default to system choosing cores
    
    // Parse command line arguments
    CoreType coreType = ANY_CORE; // Default to any core (system decides)

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-S") == 0) {
            algorithm = SEQUENTIAL;
        }
        else if (strcmp(argv[i], "-P") == 0) {
            algorithm = NAIVE_PARALLEL;
        }
        else if (strcmp(argv[i], "-V") == 0) {
            algorithm = PARALLEL_VFC;
        }
        else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            // Legacy option
            numThreads = std::atoi(argv[++i]);
            if (numThreads <= 0) {
                std::cerr << "Error: Number of threads must be positive\n";
                return EXIT_FAILURE;
            }
            std::cout << "Warning: -n is deprecated. Use -a for system-decided cores.\n";
        }
        // Sequential mode core options
        else if (strcmp(argv[i], "-p") == 0) {
            if (algorithm == SEQUENTIAL) {
                coreType = P_CORE;
            } else {
                std::cerr << "Error: -p flag is only for sequential algorithm\n";
                printUsage(argv[0]);
                return EXIT_FAILURE;
            }
        }
        else if (strcmp(argv[i], "-e") == 0) {
            if (algorithm == SEQUENTIAL) {
                coreType = E_CORE;
            } else {
                std::cerr << "Error: -e flag is only for sequential algorithm\n";
                printUsage(argv[0]);
                return EXIT_FAILURE;
            }
        }
        // Parallel mode core options
        else if (strcmp(argv[i], "-pc") == 0 && i + 1 < argc) {
            if (algorithm != SEQUENTIAL) {
                pCoreCount = std::atoi(argv[++i]);
                if (pCoreCount > 8) {
                    std::cerr << "Error: P-core count cannot exceed 8\n";
                    return EXIT_FAILURE;
                }
                useSystemCores = false;
            } else {
                std::cerr << "Error: -pc flag is only for parallel algorithms\n";
                printUsage(argv[0]);
                return EXIT_FAILURE;
            }
        }
        else if (strcmp(argv[i], "-ec") == 0 && i + 1 < argc) {
            if (algorithm != SEQUENTIAL) {
                eCoreCount = std::atoi(argv[++i]);
                if (eCoreCount > 16) {
                    std::cerr << "Error: E-core count cannot exceed 16\n";
                    return EXIT_FAILURE;
                }
                useSystemCores = false;
            } else {
                std::cerr << "Error: -ec flag is only for parallel algorithms\n";
                printUsage(argv[0]);
                return EXIT_FAILURE;
            }
        }
        else if (strcmp(argv[i], "-a") == 0) {
            if (algorithm == SEQUENTIAL) {
                coreType = ANY_CORE;
            } else {
                // For parallel, -a requires a thread count
                if (i + 1 < argc) {
                    numThreads = std::atoi(argv[++i]);
                    if (numThreads <= 0) {
                        std::cerr << "Error: Number of threads must be positive\n";
                        return EXIT_FAILURE;
                    }
                    useSystemCores = true;
                } else {
                    std::cerr << "Error: -a for parallel algorithms requires a thread count\n";
                    printUsage(argv[0]);
                    return EXIT_FAILURE;
                }
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
    switch (algorithm) {
        case SEQUENTIAL:
            std::cout << "Running sequential Louvain algorithm";
            if (coreType == P_CORE) {
                std::cout << " on P-cores\n";
            } else if (coreType == E_CORE) {
                std::cout << " on E-cores\n";
            } else {
                std::cout << " (core type not specified)\n";
            }
            louvainHierarchical(g, H, coreType);
            break;
        case NAIVE_PARALLEL:
            if (useSystemCores) {
                std::cout << "Running naive parallel Louvain algorithm with " << numThreads 
                          << " threads (system decided)\n";
                louvainParallel(g, H, numThreads);
            } else {
                std::cout << "Running naive parallel Louvain algorithm with " 
                          << pCoreCount << " P-cores and " << eCoreCount << " E-cores\n";
                louvainParallel(g, H, pCoreCount + eCoreCount, pCoreCount, eCoreCount);
            }
            break;
        case PARALLEL_VFC:
            // For now, VFC implementation still uses system-decided core allocation
            std::cout << "Running parallel Louvain algorithm with Vertex Following and Coloring using " 
                      << numThreads << " threads\n";
            louvainParallelVFC(g, H, numThreads);
            break;
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