#include "louvain_parallel_static_bl.h"
#include "louvain_seq.h"  // for fallback when numThreads == 1
#include "graph.h"
#include "hierarchy.h"
#include "core_type.h"

#include <omp.h>
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <cmath>

// Helper structure to store node index and its degree for workload balancing
struct NodeDegree {
    int index;
    double degree;
    
    bool operator<(const NodeDegree& other) const {
        return degree > other.degree; // Sort in descending order
    }
};

// Static parallel sweep with heterogeneous-aware workload distribution
static bool parallelLocalSweepStaticBL(const Graph &g, std::vector<int> &C, 
                                      int pCoreCount, int eCoreCount, double speed_ratio) {
    int n = g.n;
    double two_m = 2.0 * g.m;
    int numThreads = pCoreCount + eCoreCount;
    
    // 1) Build degree copy and community totals
    std::vector<double> k(n), comm_tot(n, 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        k[i] = g.degree[i];
    }
    for (int i = 0; i < n; i++) {
        comm_tot[C[i]] += k[i];
    }
    
    // 2) Create nodes sorted by degree for workload balancing
    std::vector<NodeDegree> sortedNodes(n);
    for (int i = 0; i < n; i++) {
        sortedNodes[i] = {i, k[i]};
    }
    std::sort(sortedNodes.begin(), sortedNodes.end());
    
    // 3) Calculate total compute capacity (P-cores = speed_ratio, E-cores = 1x)
    double totalCapacity = speed_ratio * pCoreCount + 1.0 * eCoreCount;
    std::vector<double> threadCapacity(numThreads);
    std::vector<std::vector<int>> threadNodes(numThreads);
    
    // Assign thread capacities based on core type
    for (int t = 0; t < numThreads; t++) {
        if (t < pCoreCount) {
            threadCapacity[t] = speed_ratio;  // P-cores have speed_ratio capacity
            std::cout << "Thread " << t << " is a P-core (" << speed_ratio << "x capacity)" << std::endl;
        } else {
            threadCapacity[t] = 1.0;  // E-cores have 1x capacity
            std::cout << "Thread " << t << " is an E-core (1x capacity)" << std::endl;
        }
    }
    
    // Calculate total workload (sum of all node degrees)
    double totalWorkload = 0.0;
    for (int i = 0; i < n; i++) {
        totalWorkload += k[i];
    }
    
    std::cout << "Total workload: " << totalWorkload << ", Total capacity: " << totalCapacity << std::endl;
    
    // Calculate target workload for each thread based on its capacity
    std::vector<double> targetWorkload(numThreads);
    std::vector<double> currentWorkload(numThreads, 0.0);
    
    for (int t = 0; t < numThreads; t++) {
        targetWorkload[t] = (threadCapacity[t] / totalCapacity) * totalWorkload;
        std::cout << "Thread " << t << " target workload: " << targetWorkload[t] << std::endl;
    }
    
    // 4) Distribute nodes to threads according to capacity
    for (const auto& node : sortedNodes) {
        // Find thread that is furthest from its target workload
        int bestThread = 0;
        double bestRatio = currentWorkload[0] / targetWorkload[0];
        
        for (int t = 1; t < numThreads; t++) {
            double ratio = currentWorkload[t] / targetWorkload[t];
            if (ratio < bestRatio) {
                bestRatio = ratio;
                bestThread = t;
            }
        }
        
        // Assign node to this thread
        threadNodes[bestThread].push_back(node.index);
        currentWorkload[bestThread] += node.degree;
    }
    
    // Print workload distribution
    for (int t = 0; t < numThreads; t++) {
        std::cout << "Thread " << t << " (";
        if (t < pCoreCount) std::cout << "P-core): ";
        else std::cout << "E-core): ";
        
        std::cout << "Assigned " << threadNodes[t].size() << " nodes with workload " 
                  << currentWorkload[t] << " (" 
                  << (currentWorkload[t] / targetWorkload[t] * 100) << "% of target)" 
                  << std::endl;
    }
    
    // 5) Compute best moves in parallel with static assignment
    std::vector<int> bestC(n);
    bool moved = false;
    
    #pragma omp parallel reduction(||:moved)
    {
        int tid = omp_get_thread_num();
        const std::vector<int>& myNodes = threadNodes[tid];
        
        for (int i : myNodes) {
            int curr = C[i];
            double ki = k[i];
            
            // Sum weights to each neighbor community
            std::unordered_map<int,double> w2c;
            for (auto &e : g.adj[i]) {
                w2c[C[e.neighbor]] += e.weight;
            }
            
            // Find best community
            double bestGain = 0.0;
            int bestComm = curr;
            for (auto &p : w2c) {
                double gain = p.second - (comm_tot[p.first] * ki) / two_m;
                if (gain > bestGain) {
                    bestGain = gain;
                    bestComm = p.first;
                }
            }
            bestC[i] = bestComm;
            if (bestComm != curr) moved = true;
        }
    }
    
    // 6) Apply all moves (sequentially)
    for (int i = 0; i < n; i++) {
        C[i] = bestC[i];
    }
    return moved;
}

// Aggregate graph by community vector (same as in louvain_parallel.cpp)
static Graph aggregateGraph(const Graph &g,
                            const std::vector<int> &C,
                            std::vector<int> &old2new) {
    int n = g.n;
    std::unordered_map<int,int> remap;
    old2new.resize(n);
    int cid = 0;
    // assign new IDs
    for (int i = 0; i < n; i++) {
        auto it = remap.find(C[i]);
        if (it == remap.end()) {
            remap[C[i]] = cid++;
        }
        old2new[i] = remap[C[i]];
    }
    // build new graph
    Graph g2;
    g2.n = cid;
    g2.adj.assign(cid, {});
    g2.degree.assign(cid, 0.0);

    std::vector<std::unordered_map<int,double>> W(cid);
    for (int i = 0; i < n; i++) {
        int ci = old2new[i];
        for (auto &e : g.adj[i]) {
            int cj = old2new[e.neighbor];
            W[ci][cj] += e.weight;
        }
    }
    double tot = 0.0;
    for (int i = 0; i < cid; i++) {
        for (auto &p : W[i]) {
            int j = p.first;
            double w = p.second;
            if (w > 0) {
                g2.adj[i].push_back({j, w});
                g2.degree[i] += w;
                if (i <= j) tot += w;
            }
        }
    }
    g2.m = tot;
    return g2;
}

// Top‐level parallel Louvain with heterogeneous-aware static scheduling
void louvainParallelStaticBL(const Graph &g0, Hierarchy &H, int numThreads, int pCoreCount, int eCoreCount, double speed_ratio) {
    // Fast fallback if only 1 thread requested
    if (numThreads <= 1) {
        louvainHierarchical(g0, H);
        return;
    }
    
    // Validate core counts
    if (pCoreCount + eCoreCount == 0) {
        // If no specific core counts provided, distribute threads to P and E cores
        // assuming half P cores, half E cores
        pCoreCount = numThreads / 2;
        eCoreCount = numThreads - pCoreCount;
        std::cout << "No core counts specified. Using " << pCoreCount << " P-cores and " 
                  << eCoreCount << " E-cores" << std::endl;
    } else {
        // Use provided counts, but adjust numThreads
        numThreads = pCoreCount + eCoreCount;
        std::cout << "Using " << pCoreCount << " P-cores and " << eCoreCount 
                  << " E-cores (total: " << numThreads << " threads)" << std::endl;
    }
    
    std::cout << "Speed ratio (P:E): " << speed_ratio << ":1" << std::endl;
    
    // Set thread count
    omp_set_num_threads(numThreads);
    
    // Create vector of core assignments
    std::vector<int> coreAssignments = assignParallelCores(pCoreCount, eCoreCount);
    
    if (coreAssignments.empty()) {
        // Error already reported in assignParallelCores
        return;
    }
    
    // Print the core assignments
    std::cout << "Core assignments: ";
    for (size_t i = 0; i < coreAssignments.size(); i++) {
        std::cout << coreAssignments[i];
        if (i < coreAssignments.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
    
    // Set thread affinity ONCE at the beginning
    #pragma omp parallel
    {
        int threadId = omp_get_thread_num();
        if (threadId < static_cast<int>(coreAssignments.size())) {
            setThreadAffinityToCpu(coreAssignments[threadId]);
            
            // Print confirmation from each thread
            #pragma omp critical
            {
                std::cout << "Thread " << threadId << " pinned to CPU " 
                          << coreAssignments[threadId] 
                          << " (" << (threadId < pCoreCount ? "P-core" : "E-core") << ")"
                          << std::endl;
            }
        }
    }

    // Main algorithm
    Graph g = g0;
    std::vector<int> C(g.n);
    std::iota(C.begin(), C.end(), 0);
    H.partitions.clear();

    const int MAX_LEVELS = 10;
    const int MAX_SWEEPS = 100;
    const double EPS = 1e-6;

    for (int level = 0; level < MAX_LEVELS; level++) {
        double curQ = computeModularity(g, C);
        std::cout << "Level " << level << " starting modularity: " << curQ << std::endl;

        // Local‐move sweeps with heterogeneous-aware static scheduling
        for (int sweep = 0; sweep < MAX_SWEEPS; sweep++) {
            bool moved = parallelLocalSweepStaticBL(g, C, pCoreCount, eCoreCount, speed_ratio);
            double newQ = computeModularity(g, C);
            std::cout << "  Sweep " << sweep << ": modularity = " << newQ << std::endl;
            
            if (!moved || std::fabs(newQ - curQ) < EPS) {
                break;
            }
            curQ = newQ;
        }

        // Record partition at this level
        H.partitions.push_back(C);
        std::cout << "Level " << level << " final modularity: " << curQ << std::endl;

        // Aggregate for next level
        std::vector<int> old2new;
        Graph g2 = aggregateGraph(g, C, old2new);
        if (g2.n == g.n || g2.n <= 1) {
            std::cout << "No further aggregation possible. Stopping." << std::endl;
            break;
        }
        
        std::cout << "Aggregated graph: " << g.n << " -> " << g2.n << " nodes" << std::endl;
        g = std::move(g2);
        C.assign(g.n, 0);
        std::iota(C.begin(), C.end(), 0);
    }
}