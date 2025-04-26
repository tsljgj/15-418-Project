#include "louvain_parallel.h"
#include "louvain_seq.h"  // for fallback when numThreads == 1
#include "graph.h"
#include "hierarchy.h"

#include <omp.h>
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <cmath>

// --------- Parallel sweep: one OpenMP for each node ---------
static bool parallelLocalSweep(const Graph &g, std::vector<int> &C) {
    int n = g.n;
    double two_m = 2.0 * g.m;

    // 1) Build degree copy and community totals
    std::vector<double> k(n), comm_tot(n, 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        k[i] = g.degree[i];
    }
    for (int i = 0; i < n; i++) {
        comm_tot[C[i]] += k[i];
    }

    // 2) Prepare random order for this sweep
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), std::mt19937{std::random_device{}()});

    // 3) Compute best moves in parallel
    std::vector<int> bestC(n);
    bool moved = false;
    #pragma omp parallel for schedule(dynamic,64) reduction(||:moved)
    for (int idx = 0; idx < n; idx++) {
        int i = order[idx];
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

    // 4) Apply all moves (sequentially)
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

// --------- Top‐level parallel Louvain ---------
void louvainParallel(const Graph &g0, Hierarchy &H, int numThreads) {
    // Fast fallback if only 1 thread requested
    if (numThreads <= 1) {
        louvainHierarchical(g0, H);
        return;
    }
    omp_set_num_threads(numThreads);

    Graph g = g0;
    std::vector<int> C(g.n);
    std::iota(C.begin(), C.end(), 0);
    H.partitions.clear();

    const int MAX_LEVELS = 10;
    const int MAX_SWEEPS = 100;
    const double EPS    = 1e-6;

    for (int level = 0; level < MAX_LEVELS; level++) {
        double curQ = computeModularity(g, C);

        // Local‐move sweeps
        for (int sweep = 0; sweep < MAX_SWEEPS; sweep++) {
            bool moved = parallelLocalSweep(g, C);
            double newQ = computeModularity(g, C);
            if (!moved || std::fabs(newQ - curQ) < EPS) {
                break;
            }
            curQ = newQ;
        }

        // Record partition at this level
        H.partitions.push_back(C);

        // Aggregate for next level
        std::vector<int> old2new;
        Graph g2 = aggregateGraph(g, C, old2new);
        if (g2.n == g.n || g2.n <= 1) break;
        g = std::move(g2);
        C.assign(g.n, 0);
        std::iota(C.begin(), C.end(), 0);
    }
}