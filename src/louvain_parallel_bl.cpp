#include "louvain_parallel_bl.h"
#include "louvain_seq.h"
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
#include <atomic>

// --------- P-core work: compute-intensive tasks ---------
static void computeModularityGains(const Graph &g, const std::vector<int> &C,
                                   const std::vector<double> &k, 
                                   const std::vector<double> &comm_tot,
                                   int start, int end, double two_m,
                                   std::vector<int> &bestC,
                                   std::atomic<bool> &moved) {
    for (int idx = start; idx < end; idx++) {
        int i = idx;
        int curr = C[i];
        double ki = k[i];

        // Sum weights to each neighbor community
        std::unordered_map<int,double> w2c;
        for (auto &e : g.adj[i]) {
            w2c[C[e.neighbor]] += e.weight;
        }

        // Find best community (computationally intensive)
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

// --------- E-core work: memory-intensive tasks ---------
static void prepareDataStructures(const Graph &g, const std::vector<int> &C,
                                  int start, int end,
                                  std::vector<double> &k_local) {
    for (int i = start; i < end; i++) {
        k_local[i] = g.degree[i];
    }
}

// --------- Parallel sweep with P/E core specialization ---------
static bool parallelLocalSweep(const Graph &g, std::vector<int> &C, 
                               int pCoreCount, int eCoreCount) {
    int n = g.n;
    double two_m = 2.0 * g.m;

    // Data structures
    std::vector<double> k(n);
    std::vector<double> comm_tot(n, 0.0);

    // Phase 1: E-cores handle data preparation
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int total_threads = pCoreCount + eCoreCount;
        bool is_e_core = (tid >= pCoreCount);

        if (is_e_core) {
            // E-cores prepare degree array
            int chunk_size = n / eCoreCount;
            int start = (tid - pCoreCount) * chunk_size;
            int end = (tid == total_threads - 1) ? n : start + chunk_size;
            
            prepareDataStructures(g, C, start, end, k);
        }
    }
    
    #pragma omp barrier

    // Phase 2: E-cores update community totals with OpenMP atomic
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int total_threads = pCoreCount + eCoreCount;
        bool is_e_core = (tid >= pCoreCount);

        if (is_e_core) {
            int chunk_size = n / eCoreCount;
            int start = (tid - pCoreCount) * chunk_size;
            int end = (tid == total_threads - 1) ? n : start + chunk_size;
            
            for (int i = start; i < end; i++) {
                #pragma omp atomic
                comm_tot[C[i]] += k[i];
            }
        }
    }

    #pragma omp barrier

    // Phase 3: P-cores compute best moves
    std::vector<int> bestC(n);
    std::atomic<bool> moved(false);

    // Categorize nodes by computational intensity
    std::vector<int> dense_nodes, sparse_nodes;
    for (int i = 0; i < n; i++) {
        if (g.adj[i].size() > 10) { // Threshold for "dense"
            dense_nodes.push_back(i);
        } else {
            sparse_nodes.push_back(i);
        }
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        bool is_p_core = (tid < pCoreCount);

        if (is_p_core) {
            // P-cores handle dense nodes (more compute-intensive)
            int chunk_size = dense_nodes.size() / pCoreCount;
            int start = tid * chunk_size;
            int end = (tid == pCoreCount - 1) ? dense_nodes.size() : start + chunk_size;
            
            for (int idx = start; idx < end; idx++) {
                int i = dense_nodes[idx];
                computeModularityGains(g, C, k, comm_tot, i, i+1, two_m, bestC, moved);
            }
        } else {
            // E-cores handle sparse nodes (less compute-intensive)
            int e_core_id = tid - pCoreCount;
            int chunk_size = sparse_nodes.size() / eCoreCount;
            int start = e_core_id * chunk_size;
            int end = (e_core_id == eCoreCount - 1) ? sparse_nodes.size() : start + chunk_size;
            
            for (int idx = start; idx < end; idx++) {
                int i = sparse_nodes[idx];
                computeModularityGains(g, C, k, comm_tot, i, i+1, two_m, bestC, moved);
            }
        }
    }

    // Phase 4: Apply all moves (E-cores handle this memory operation)
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        bool is_e_core = (tid >= pCoreCount);
        int total_threads = pCoreCount + eCoreCount;

        if (is_e_core) {
            int chunk_size = n / eCoreCount;
            int start = (tid - pCoreCount) * chunk_size;
            int end = (tid == total_threads - 1) ? n : start + chunk_size;
            
            for (int i = start; i < end; i++) {
                C[i] = bestC[i];
            }
        }
    }

    return moved.load();
}

// --------- Aggregate graph with P/E core specialization ---------
static Graph aggregateGraph(const Graph &g, const std::vector<int> &C,
                            std::vector<int> &old2new, int pCoreCount, int eCoreCount) {
    int n = g.n;
    std::unordered_map<int,int> remap;
    old2new.resize(n);
    int cid = 0;

    // P-cores handle community remapping (more complex logic)
    for (int i = 0; i < n; i++) {
        auto it = remap.find(C[i]);
        if (it == remap.end()) {
            remap[C[i]] = cid++;
        }
        old2new[i] = remap[C[i]];
    }

    Graph g2;
    g2.n = cid;
    g2.adj.assign(cid, {});
    g2.degree.assign(cid, 0.0);

    // Use thread-local storage for parallel edge aggregation
    std::vector<std::vector<std::unordered_map<int,double>>> thread_local_W;
    int total_threads = pCoreCount + eCoreCount;
    thread_local_W.resize(total_threads);
    for (int t = 0; t < total_threads; t++) {
        thread_local_W[t].resize(cid);
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        bool is_p_core = (tid < pCoreCount);

        if (is_p_core) {
            // P-cores process dense regions
            #pragma omp for schedule(dynamic, 64)
            for (int i = 0; i < n; i++) {
                if (g.adj[i].size() > 10) { // Dense nodes
                    int ci = old2new[i];
                    for (auto &e : g.adj[i]) {
                        int cj = old2new[e.neighbor];
                        thread_local_W[tid][ci][cj] += e.weight;
                    }
                }
            }
        } else {
            // E-cores process sparse regions
            #pragma omp for schedule(dynamic, 128)
            for (int i = 0; i < n; i++) {
                if (g.adj[i].size() <= 10) { // Sparse nodes
                    int ci = old2new[i];
                    for (auto &e : g.adj[i]) {
                        int cj = old2new[e.neighbor];
                        thread_local_W[tid][ci][cj] += e.weight;
                    }
                }
            }
        }
    }

    // Merge thread-local results (E-cores handle this memory operation)
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        bool is_e_core = (tid >= pCoreCount);

        if (is_e_core) {
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < cid; i++) {
                std::unordered_map<int,double> merged;
                for (int t = 0; t < total_threads; t++) {
                    for (auto &p : thread_local_W[t][i]) {
                        merged[p.first] += p.second;
                    }
                }
                
                for (auto &p : merged) {
                    int j = p.first;
                    double w = p.second;
                    if (w > 0) {
                        g2.adj[i].push_back({j, w});
                        g2.degree[i] += w;
                    }
                }
            }
        }
    }

    // Final edge count (single threaded for correctness)
    double tot = 0.0;
    for (int i = 0; i < cid; i++) {
        for (auto &e : g2.adj[i]) {
            if (i <= e.neighbor) tot += e.weight;
        }
    }
    g2.m = tot;
    
    return g2;
}

// --------- Top‐level parallel Louvain with P/E optimization ---------
void louvainParallelBL(const Graph &g0, Hierarchy &H, int numThreads, 
                     int pCoreCount, int eCoreCount) {
    // Fast fallback if only 1 thread requested
    if (numThreads <= 1) {
        louvainHierarchical(g0, H);
        return;
    }
    
    // Set up core allocation
    std::vector<int> coreAssignments;
    bool useSpecificCores = (pCoreCount > 0 || eCoreCount > 0);
    
    if (useSpecificCores) {
        coreAssignments = assignParallelCores(pCoreCount, eCoreCount);
        
        if (coreAssignments.empty()) {
            return;
        }
        
        numThreads = pCoreCount + eCoreCount;
        
        std::cout << "Running with " << pCoreCount << " P-cores and " 
                  << eCoreCount << " E-cores\n";
    } else {
        std::cout << "Using " << numThreads << " threads with system-decided core affinity\n";
    }
    
    omp_set_num_threads(numThreads);
    
    // Set thread affinity
    if (useSpecificCores) {
        #pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            if (threadId < static_cast<int>(coreAssignments.size())) {
                setThreadAffinityToCpu(coreAssignments[threadId]);
            }
        }
    }

    Graph g = g0;
    std::vector<int> C(g.n);
    std::iota(C.begin(), C.end(), 0);
    H.partitions.clear();

    const int MAX_LEVELS = 10;
    const int MAX_SWEEPS = 100;
    const double EPS = 1e-6;

    for (int level = 0; level < MAX_LEVELS; level++) {
        double curQ = computeModularity(g, C);

        // Local‐move sweeps with P/E optimization
        for (int sweep = 0; sweep < MAX_SWEEPS; sweep++) {
            bool moved = useSpecificCores ? 
                         parallelLocalSweep(g, C, pCoreCount, eCoreCount) :
                         parallelLocalSweep(g, C, numThreads, 0);
            double newQ = computeModularity(g, C);
            if (!moved || std::fabs(newQ - curQ) < EPS) {
                break;
            }
            curQ = newQ;
        }

        H.partitions.push_back(C);

        // Aggregate for next level
        std::vector<int> old2new;
        Graph g2 = useSpecificCores ?
                   aggregateGraph(g, C, old2new, pCoreCount, eCoreCount) :
                   aggregateGraph(g, C, old2new, numThreads, 0);
        
        if (g2.n == g.n || g2.n <= 1) break;
        
        g = std::move(g2);
        C.assign(g.n, 0);
        std::iota(C.begin(), C.end(), 0);
    }
}