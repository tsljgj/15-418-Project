// louvain_parallel_bl.cpp
// Big-LITTLE aware parallel Louvain implementation

#include "louvain_parallel_bl.h"
#include "louvain_seq.h"
#include "graph.h"
#include "hierarchy.h"
#include "core_type.h"

#include <omp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <atomic>

// --- graph aggregation (same as in naive version) ---
static Graph aggregateGraph(const Graph &g,
                            const std::vector<int> &C,
                            std::vector<int> &old2new) {
    int n = g.n;
    std::unordered_map<int,int> remap;
    old2new.resize(n);
    int cid = 0;
    for (int i = 0; i < n; ++i) {
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

    std::vector<std::unordered_map<int,double>> W(cid);
    for (int i = 0; i < n; ++i) {
        int ci = old2new[i];
        for (auto &e : g.adj[i]) {
            int cj = old2new[e.neighbor];
            W[ci][cj] += e.weight;
        }
    }
    double tot = 0.0;
    for (int i = 0; i < cid; ++i) {
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

void louvainParallelBL(const Graph &g0, Hierarchy &H,
                       int numThreads, int pCoreCount, int eCoreCount) {
    if (numThreads <= 1) {
        // Fallback to sequential hierarchical Louvain
        louvainHierarchical(g0, H);
        return;
    }

    bool useSpecificCores = (pCoreCount > 0 || eCoreCount > 0);
    std::vector<int> coreAssignments;

    // Assign P/E cores if requested
    if (useSpecificCores) {
        coreAssignments = assignParallelCores(pCoreCount, eCoreCount);
        if (coreAssignments.empty()) return;
        numThreads = (int)coreAssignments.size();
        // std::cout << "Core assignments: ";
        // for (int id : coreAssignments) std::cout << id << " ";
        // std::cout << std::endl;
    } else {
        std::cout << "Using " << numThreads << " threads (auto affinity)\n";
    }

    omp_set_num_threads(numThreads);

    // Pin threads once before sweeping
    if (useSpecificCores) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            if (tid < (int)coreAssignments.size()) {
                setThreadAffinityToCpu(coreAssignments[tid]);
                // #pragma omp critical
                // std::cout << "Thread " << tid
                //           << " pinned to CPU " << coreAssignments[tid] << "\n";
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

    for (int level = 0; level < MAX_LEVELS; ++level) {
        double curQ = computeModularity(g, C);

        for (int sweep = 0; sweep < MAX_SWEEPS; ++sweep) {
            int n = g.n;
            double two_m = 2.0 * g.m;

            // 1) Build degree copy and community totals
            std::vector<double> k = g.degree;
            std::vector<double> comm_tot(n, 0.0);
            for (int i = 0; i < n; ++i)
                comm_tot[C[i]] += k[i];

            // 2) Partition nodes by "intensity" (degree-based)
            double avg_deg = std::accumulate(k.begin(), k.end(), 0.0) / n;
            std::vector<int> heavyNodes, lightNodes;
            heavyNodes.reserve(n);
            lightNodes.reserve(n);
            for (int i = 0; i < n; ++i) {
                if (k[i] > avg_deg) heavyNodes.push_back(i);
                else               lightNodes.push_back(i);
            }

            // 3) Two-queue work-stealing
            std::atomic<int> hPtr(0), lPtr(0);
            std::atomic<bool> moved(false);
            std::vector<int> bestC(n);

            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                bool isP = useSpecificCores ? (tid < pCoreCount) : true;

                while (true) {
                    int idx;

                    // P-threads drain heavy first
                    if (isP) {
                        idx = hPtr.fetch_add(1);
                        if (idx < (int)heavyNodes.size()) {
                            int i = heavyNodes[idx];
                            int curr = C[i];
                            double ki = k[i];
                            std::unordered_map<int,double> w2c;
                            for (auto &e : g.adj[i])
                                w2c[C[e.neighbor]] += e.weight;
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
                            if (bestComm != curr) moved.store(true);
                            continue;
                        }
                        // heavy empty → fall through to light
                    }

                    // Drain light nodes
                    idx = lPtr.fetch_add(1);
                    if (idx < (int)lightNodes.size()) {
                        int i = lightNodes[idx];
                        int curr = C[i];
                        double ki = k[i];
                        std::unordered_map<int,double> w2c;
                        for (auto &e : g.adj[i])
                            w2c[C[e.neighbor]] += e.weight;
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
                        if (bestComm != curr) moved.store(true);
                        continue;
                    }

                    // E-threads can steal heavy when done with light
                    if (!isP) {
                        idx = hPtr.fetch_add(1);
                        if (idx < (int)heavyNodes.size()) {
                            int i = heavyNodes[idx];
                            int curr = C[i];
                            double ki = k[i];
                            std::unordered_map<int,double> w2c;
                            for (auto &e : g.adj[i])
                                w2c[C[e.neighbor]] += e.weight;
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
                            if (bestComm != curr) moved.store(true);
                            continue;
                        }
                    }

                    // No more work
                    break;
                }
            } // end parallel

            // 4) Apply all moves
            for (int i = 0; i < n; ++i)
                C[i] = bestC[i];

            double newQ = computeModularity(g, C);
            if (!moved.load() || std::fabs(newQ - curQ) < EPS)
                break;
            curQ = newQ;
        }

        // Record partition and aggregate
        H.partitions.push_back(C);
        std::vector<int> old2new;
        Graph g2 = aggregateGraph(g, C, old2new);
        if (g2.n == g.n || g2.n <= 1) break;
        g = std::move(g2);
    }
}
