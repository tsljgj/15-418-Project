// louvain_static.cpp

#include "louvain_parallel_static.h"
#include "louvain_seq.h"       // for louvainHierarchical fallback
#include "graph.h"
#include "hierarchy.h"
#include "core_type.h"         // assignParallelCores, setThreadAffinityToCpu

#include <omp.h>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cmath>

// --------- Degree‐Balanced Static Local Sweep ---------
static bool parallelLocalSweepStatic(const Graph &g, std::vector<int> &C) {
    int n       = g.n;
    double two_m = 2.0 * g.m;

    // 1) Build degree copy and community totals (serial)
    std::vector<double> k(n), comm_tot(n, 0.0);
    for (int i = 0; i < n; i++) {
        k[i] = g.degree[i];
    }
    for (int i = 0; i < n; i++) {
        comm_tot[C[i]] += k[i];
    }

    // 2) Partition nodes into T bins to balance sum of k[i]
    int T = omp_get_max_threads();
    std::vector<std::vector<int>> bins(T);
    std::vector<double> loads(T, 0.0);

    // Sort nodes by descending degree
    std::vector<int> nodes(n);
    std::iota(nodes.begin(), nodes.end(), 0);
    std::sort(nodes.begin(), nodes.end(),
              [&](int a, int b){ return k[a] > k[b]; });

    // Greedily assign each node to the lightest‐loaded bin
    for (int i : nodes) {
        int t = std::min_element(loads.begin(), loads.end()) - loads.begin();
        bins[t].push_back(i);
        loads[t] += k[i];
    }

    // 3) Each thread processes its own bin
    std::vector<int> bestC(n);
    bool moved = false;
    #pragma omp parallel reduction(||:moved)
    {
        int tid = omp_get_thread_num();
        for (int i : bins[tid]) {
            int curr = C[i];
            double ki = k[i];

            // sum weights to each neighboring community
            std::unordered_map<int,double> w2c;
            for (auto &e : g.adj[i]) {
                w2c[C[e.neighbor]] += e.weight;
            }

            // find best community move
            double bestGain = 0.0;
            int    bestComm = curr;
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

    // 4) Apply all moves
    for (int i = 0; i < n; i++) {
        C[i] = bestC[i];
    }

    return moved;
}

// --------- Aggregate graph by community vector ---------
static Graph aggregateGraph(const Graph &g,
                            const std::vector<int> &C,
                            std::vector<int> &old2new) {
    int n = g.n;
    std::unordered_map<int,int> remap;
    old2new.resize(n);
    int cid = 0;
    for (int i = 0; i < n; i++) {
        auto it = remap.find(C[i]);
        if (it == remap.end()) remap[C[i]] = cid++;
        old2new[i] = remap[C[i]];
    }

    Graph g2;
    g2.n      = cid;
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
            int j    = p.first;
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

// --------- Top‐level static‐schedule Louvain ---------
void louvainParallelStatic(const Graph &g0,
                           Hierarchy &H,
                           int numThreads,
                           int pCoreCount,
                           int eCoreCount) {
    // Fast fallback if only 1 thread requested
    if (numThreads <= 1) {
        louvainHierarchical(g0, H);
        return;
    }
    
    // Create vector of core assignments
    std::vector<int> coreAssignments;
    bool useSpecificCores = (pCoreCount > 0 || eCoreCount > 0);
    
    // If specific core counts are provided, set up thread affinity
    if (useSpecificCores) {
        coreAssignments = assignParallelCores(pCoreCount, eCoreCount);
        
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
        
        // The total number of threads will be the sum of the core counts
        numThreads = pCoreCount + eCoreCount;
    } else {
        std::cout << "Using " << numThreads << " threads with system-decided core affinity\n";
    }
    
    omp_set_num_threads(numThreads);
    
    // Set thread affinity ONCE at the beginning
    if (useSpecificCores) {
        #pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            if (threadId < static_cast<int>(coreAssignments.size())) {
                setThreadAffinityToCpu(coreAssignments[threadId]);
                
                // Print confirmation from each thread
                #pragma omp critical
                {
                    std::cout << "Thread " << threadId << " pinned to CPU " 
                              << coreAssignments[threadId] << std::endl;
                }
            }
        }
        
    }

    // 5) Initialize clustering
    Graph g = g0;
    std::vector<int> C(g.n);
    std::iota(C.begin(), C.end(), 0);
    H.partitions.clear();

    const int MAX_LEVELS = 10;
    const int MAX_SWEEPS = 100;
    const double EPS    = 1e-6;

    // 6) Hierarchical Louvain
    for (int lvl = 0; lvl < MAX_LEVELS; lvl++) {
        double curQ = computeModularity(g, C);

        // Local‐move sweeps
        for (int sw = 0; sw < MAX_SWEEPS; sw++) {
            bool moved = parallelLocalSweepStatic(g, C);
            double newQ = computeModularity(g, C);
            if (!moved || std::fabs(newQ - curQ) < EPS) break;
            curQ = newQ;
        }

        // record partition
        H.partitions.push_back(C);

        // coarsen
        std::vector<int> old2new;
        Graph g2 = aggregateGraph(g, C, old2new);
        if (g2.n == g.n || g2.n <= 1) break;
        g = std::move(g2);
        C.assign(g.n, 0);
        std::iota(C.begin(), C.end(), 0);
    }
}
