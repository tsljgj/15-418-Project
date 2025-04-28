// louvain_naive_bl.cpp
#include "louvain_parallel_bl.h"
#include "louvain_seq.h"      // for fallback
#include "graph.h"
#include "hierarchy.h"
#include "core_type.h"

#include <omp.h>
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <cmath>

// inside louvain_naive_bl.cpp, replace parallelLocalSweepBL with:

static bool parallelLocalSweepBL(
    const Graph &g,
    std::vector<int> &C,
    int pCoreCount,
    int eCoreCount
) {
    int n = g.n;
    double two_m = 2.0 * g.m;

    // 1) Copy degrees + community totals
    std::vector<double> k(n), comm_tot(n, 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        k[i] = g.degree[i];
    }
    for (int i = 0; i < n; i++) {
        comm_tot[C[i]] += k[i];
    }

    // 2) Build a list of all nodes sorted by descending degree
    std::vector<int> sorted(n);
    std::iota(sorted.begin(), sorted.end(), 0);
    std::sort(sorted.begin(), sorted.end(),
              [&](int a, int b){ return k[a] > k[b]; });

    // 3) Compute how many heavy nodes go to P-cores
    int totalThreads = pCoreCount + eCoreCount;
    int heavyCount  = int((double)n * pCoreCount / totalThreads);

    // 4) Build per-thread static lists
    std::vector<std::vector<int>> threadNodes(totalThreads);
    // 4a) Top heavyCount → P-cores (threads 0…pCoreCount-1)
    for (int idx = 0; idx < heavyCount; ++idx) {
        int tid = idx % pCoreCount;
        threadNodes[tid].push_back(sorted[idx]);
    }
    // 4b) Remaining → E-cores (threads pCoreCount…p+e-1)
    for (int idx = heavyCount; idx < n; ++idx) {
        int tid = pCoreCount + ((idx - heavyCount) % eCoreCount);
        threadNodes[tid].push_back(sorted[idx]);
    }

    // 5) Each thread processes its own nodes
    std::vector<int> bestC(n);
    bool moved = false;
    #pragma omp parallel reduction(||:moved)
    {
        int tid = omp_get_thread_num();
        auto &myList = threadNodes[tid];

        // thread-local temp map
        std::unordered_map<int,double> w2c;
        for (int i : myList) {
            int curr = C[i];
            double ki = k[i];
            w2c.clear();

            // accumulate weights to neighbor communities
            for (auto &e : g.adj[i]) {
                w2c[C[e.neighbor]] += e.weight;
            }

            // find best move
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

    // 6) Apply all moves
    for (int i = 0; i < n; i++) {
        C[i] = bestC[i];
    }
    return moved;
}


// reuse the same aggregation as naive
static Graph aggregateGraph(
    const Graph &g,
    const std::vector<int> &C,
    std::vector<int> &old2new
) {
    int n = g.n;
    std::unordered_map<int,int> remap;
    old2new.resize(n);
    int cid = 0;
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

// --------- Top-level big-little Louvain ---------
void louvainParallelBL(
    const Graph &g0,
    Hierarchy &H,
    int numThreads,
    int pCoreCount,
    int eCoreCount
) {
    if (numThreads <= 1) {
        louvainHierarchical(g0, H);
        return;
    }

    // assign specific cores
    auto coreAssignments = assignParallelCores(pCoreCount, eCoreCount);
    if (coreAssignments.empty()) {
        return; // error already printed
    }
    numThreads = pCoreCount + eCoreCount;
    omp_set_num_threads(numThreads);

    // bind each thread once
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        if (tid < (int)coreAssignments.size()) {
            setThreadAffinityToCpu(coreAssignments[tid]);
        }
    }

    // iterative louvain
    Graph g = g0;
    std::vector<int> C(g.n);
    std::iota(C.begin(), C.end(), 0);
    H.partitions.clear();

    const int MAX_LEVELS = 10;
    const int MAX_SWEEPS = 100;
    const double EPS    = 1e-6;

    for (int level = 0; level < MAX_LEVELS; level++) {
        double curQ = computeModularity(g, C);

        for (int sweep = 0; sweep < MAX_SWEEPS; sweep++) {
            bool moved = parallelLocalSweepBL(g, C, pCoreCount, eCoreCount);
            double newQ = computeModularity(g, C);
            if (!moved || std::fabs(newQ - curQ) < EPS) {
                break;
            }
            curQ = newQ;
        }

        H.partitions.push_back(C);
        std::vector<int> old2new;
        Graph g2 = aggregateGraph(g, C, old2new);
        if (g2.n == g.n || g2.n <= 1) break;
        g = std::move(g2);
        C.assign(g.n, 0);
        std::iota(C.begin(), C.end(), 0);
    }
}
