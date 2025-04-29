// louvain_parallel_static_bl.cpp

#include "louvain_parallel_static_bl.h"
#include "louvain_seq.h"     // for louvainHierarchical fallback
#include "graph.h"
#include "hierarchy.h"

#include <omp.h>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cmath>

// --------- Heterogeneity‐aware Static Local Sweep ---------
static bool parallelLocalSweepBL(const Graph &g,
                                 std::vector<int> &C,
                                 int pCoreCount,
                                 int eCoreCount)
{
    int n       = g.n;
    double two_m = 2.0 * g.m;

    // 1) Copy degrees and accumulate community totals
    std::vector<double> k(n), comm_tot(n, 0.0);
    for (int i = 0; i < n; i++) {
        k[i] = g.degree[i];
    }
    for (int i = 0; i < n; i++) {
        comm_tot[C[i]] += k[i];
    }

    // 2) Build speed‐weighted bins
    int T = pCoreCount + eCoreCount;
    std::vector<std::vector<int>> bins(T);
    std::vector<double> loads(T, 0.0);

    // Relative core speeds
    const double pSpeed = 3.0, eSpeed = 1.0;
    std::vector<double> coreSpeed(T);
    for (int t = 0; t < T; t++) {
        coreSpeed[t] = (t < pCoreCount ? pSpeed : eSpeed);
    }

    // 2a) Sort nodes by descending degree
    std::vector<int> nodes(n);
    std::iota(nodes.begin(), nodes.end(), 0);
    std::sort(nodes.begin(), nodes.end(),
              [&](int a, int b){ return k[a] > k[b]; });

    // 2b) Greedy binning by estimated time = degree / speed
    for (int i : nodes) {
        int bestT = 0;
        double bestLoad = loads[0];
        for (int t = 1; t < T; t++) {
            if (loads[t] < bestLoad) {
                bestLoad = loads[t];
                bestT = t;
            }
        }
        bins[bestT].push_back(i);
        loads[bestT] += k[i] / coreSpeed[bestT];
    }

    // 3) Parallel local sweep
    std::vector<int> bestC(n);
    bool moved = false;
    #pragma omp parallel reduction(||:moved)
    {
        int tid = omp_get_thread_num();
        for (int i : bins[tid]) {
            int curr = C[i];
            double ki = k[i];

            std::unordered_map<int,double> w2c;
            for (auto &e : g.adj[i]) {
                w2c[C[e.neighbor]] += e.weight;
            }

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

    // 4) Apply moves
    for (int i = 0; i < n; i++) {
        C[i] = bestC[i];
    }
    return moved;
}

// --------- Aggregate graph by community vector ---------
static Graph aggregateGraph(const Graph &g,
                            const std::vector<int> &C,
                            std::vector<int> &old2new)
{
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

// --------- Entry point matching the header ---------
void louvainParallelStaticBL(const Graph &g0,
                             Hierarchy &H,
                             int numThreads,
                             int pCoreCount,
                             int eCoreCount)
{
    if (numThreads <= 1) {
        louvainHierarchical(g0, H);
        return;
    }

    int T = (pCoreCount + eCoreCount > 0
             ? pCoreCount + eCoreCount
             : numThreads);
    omp_set_num_threads(T);

    Graph g = g0;
    std::vector<int> C(g.n);
    std::iota(C.begin(), C.end(), 0);
    H.partitions.clear();

    const int MAX_LEVELS = 10;
    const int MAX_SWEEPS = 100;
    const double EPS    = 1e-6;

    for (int lvl = 0; lvl < MAX_LEVELS; lvl++) {
        double curQ = computeModularity(g, C);

        for (int sw = 0; sw < MAX_SWEEPS; sw++) {
            bool moved = parallelLocalSweepBL(g, C, pCoreCount, eCoreCount);
            double newQ = computeModularity(g, C);
            if (!moved || std::fabs(newQ - curQ) < EPS) break;
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
