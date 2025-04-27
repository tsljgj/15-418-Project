#include "louvain_parallel_vfc.h"
#include "graph.h"

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <omp.h>

// -----------------------------------------------------------------------------
// 1) high-quality, balanced distance-1 coloring
// -----------------------------------------------------------------------------
std::vector<std::vector<int>> createColorSets(const Graph &g) {
    int n = g.n;
    std::vector<int> color(n, -1);
    std::vector<std::vector<int>> sets;
    std::vector<int> forbidden;
    forbidden.reserve(64);

    // Simple greedy 
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int v = 0; v < n; ++v) {
        forbidden.clear();
        for (auto &e : g.adj[v]) {
            int c = color[e.neighbor];
            if (c >= 0) forbidden.push_back(c);
        }
        std::sort(forbidden.begin(), forbidden.end());
        int c = 0;
        for (int fc : forbidden) {
            if (fc == c) ++c;
            else break;
        }
        color[v] = c;
    }
    // gather into sets
    for (int v = 0; v < n; ++v) {
        int c = color[v];
        if (c >= (int)sets.size()) sets.resize(c + 1);
        sets[c].push_back(v);
    }
    // drop empties
    sets.erase(
      std::remove_if(sets.begin(), sets.end(),
                     [](const std::vector<int> &S){ return S.empty(); }),
      sets.end());
    std::cerr << "[Color] #sets=" << sets.size() << "\n";
    return sets;
}

// -----------------------------------------------------------------------------
// 2) Vertex-following preprocessing with chain‐collapse
// -----------------------------------------------------------------------------
Graph applyVertexFollowing(const Graph &g, std::vector<int> &mapping) {
    int n = g.n;
    mapping.resize(n);
    std::iota(mapping.begin(), mapping.end(), 0);

    std::vector<bool> is_single(n,false);
    std::vector<int> nbr(n,-1);

    #pragma omp parallel for schedule(dynamic,1024)
    for (int i = 0; i < n; ++i) {
        if (g.adj[i].size() == 1) {
            is_single[i] = true;
            nbr[i] = g.adj[i][0].neighbor;
        }
    }
    // if no singles, done
    if (!std::any_of(is_single.begin(), is_single.end(),
                     [](bool b){ return b; }))
        return g;

    // assign new IDs to non‐singles
    std::vector<int> old2new(n,-1);
    int new_n = 0;
    for (int i = 0; i < n; ++i)
        if (!is_single[i])
            old2new[i] = new_n++;

    // final mapping: collapse entire chains
    for (int i = 0; i < n; ++i) {
        if (is_single[i]) {
            int t = nbr[i];
            while (is_single[t]) t = nbr[t];
            mapping[i] = old2new[t];
        } else {
            mapping[i] = old2new[i];
        }
    }

    // build compacted graph
    Graph ng;
    ng.n = new_n;
    ng.adj.assign(new_n,{});
    ng.degree.assign(new_n,0.0);
    ng.m = 0.0;

    int T = omp_get_max_threads();
    std::vector<std::vector<std::unordered_map<int,double>>> localE(
      T, std::vector<std::unordered_map<int,double>>(new_n));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic,1024)
        for (int i = 0; i < n; ++i) {
            if (is_single[i]) continue;
            int ni = old2new[i];
            for (auto &e : g.adj[i]) {
                int nj = mapping[e.neighbor];
                if (ni <= nj)
                    localE[tid][ni][nj] += e.weight;
                else
                    localE[tid][nj][ni] += e.weight;
            }
        }
    }
    // merge
    for (int i = 0; i < new_n; ++i) {
        std::unordered_map<int,double> M;
        for (int t = 0; t < T; ++t)
            for (auto &kv : localE[t][i])
                M[kv.first] += kv.second;
        for (auto &kv : M) {
            int j = kv.first; double w = kv.second;
            ng.adj[i].push_back({j,w});
            if (i==j) {
                ng.degree[i]+=2*w;
                ng.m += w;
            } else {
                ng.degree[i]+=w;
                ng.m += w;
            }
        }
    }
    std::cerr << "[VF] collapsed "<< n <<"→"<< new_n <<"\n";
    return ng;
}

// -----------------------------------------------------------------------------
// 3) One local‐move phase (min‐label + thresholding)
// -----------------------------------------------------------------------------
bool localMovePhaseWithIsolateSets(
    Graph &g,
    std::vector<int> &C,
    const std::vector<std::vector<int>> &sets,
    double threshold)
{
    double inv2m = 1.0 / (2.0 * g.m);
    double Q0 = computeModularity(g,C);
    bool anyMoved = false;

    for (auto &set : sets) {
        // community‐degree totals
        std::vector<double> commTot(g.n,0.0);
        for (int v = 0; v < g.n; ++v)
            commTot[C[v]] += g.degree[v];

        bool moved = false;
        #pragma omp parallel for schedule(dynamic,1024) reduction(|:moved)
        for (size_t k = 0; k < set.size(); ++k) {
            int v = set[k];
            int curr = C[v];
            double bestGain = 0.0;
            int bestC = curr;

            // neighbor‐weights
            std::unordered_map<int,double> W;
            for (auto &e : g.adj[v])
                W[C[e.neighbor]] += e.weight;

            for (auto &kv : W) {
                int c = kv.first;
                if (c==curr) continue;
                double gain = kv.second
                              - g.degree[v]*commTot[c]*inv2m;
                if (gain > bestGain
                   || (gain==bestGain && c<bestC)) {
                    bestGain = gain;
                    bestC = c;
                }
            }
            if (bestC!=curr) {
                C[v]=bestC;
                moved = true;
            }
        }
        if (moved) anyMoved = true;
    }
    double Q1 = computeModularity(g,C);
    return anyMoved && (Q1-Q0) > threshold * std::abs(Q0);
}

// -----------------------------------------------------------------------------
// 4) Fast parallel graph aggregation
// -----------------------------------------------------------------------------
Graph aggregateGraphFast(
    const Graph &g,
    const std::vector<int> &C,
    std::unordered_map<int,int> &c2n)
{
    c2n.clear();
    int cnt=0;
    for (int v=0; v<g.n; ++v) {
        auto it = c2n.find(C[v]);
        if (it==c2n.end())
            c2n[C[v]] = cnt++;
    }

    Graph ng;
    ng.n = cnt;
    ng.adj.assign(cnt,{});
    ng.degree.assign(cnt,0.0);
    ng.m = 0.0;

    int T = omp_get_max_threads();
    std::vector<std::vector<std::unordered_map<int,double>>> localE(
      T, std::vector<std::unordered_map<int,double>>(cnt));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic,1024)
        for (int v=0; v<g.n; ++v) {
            int cv = c2n[C[v]];
            for (auto &e: g.adj[v]) {
                int cu = c2n[C[e.neighbor]];
                if (cv <= cu)
                    localE[tid][cv][cu] += e.weight;
                else
                    localE[tid][cu][cv] += e.weight;
            }
        }
    }
    for (int i=0; i<cnt; ++i) {
        std::unordered_map<int,double> M;
        for (int t=0; t<T; ++t)
            for (auto &kv: localE[t][i])
                M[kv.first] += kv.second;
        for (auto &kv: M) {
            int j=kv.first; double w=kv.second;
            ng.adj[i].push_back({j,w});
            if (i==j) {
                ng.degree[i]+=2*w;
                ng.m += w;
            } else {
                ng.degree[i]+=w;
                ng.m += w;
            }
        }
    }
    return ng;
}

// -----------------------------------------------------------------------------
// Main driver: dynamic thresholds + conditional recoloring
// -----------------------------------------------------------------------------
void louvainParallelVFC(
    const Graph &input_g,
    Hierarchy      &H,
    int              numThreads)
{
    omp_set_num_threads(numThreads);
    std::cerr<<"Parallel VFC Louvain ("<<numThreads<<" threads)\n";

    // 1) Vertex‐following
    std::vector<int> vf_map;
    Graph g = applyVertexFollowing(input_g, vf_map);

    // map back to original
    std::vector<int> hmap(input_g.n);
    std::iota(hmap.begin(), hmap.end(), 0);

    const double TH_COLORED = 1e-2;
    const double TH_FINE    = 1e-6;
    const int    MIN_RECOLOR_VERTS = 100000;

    double lastPhaseGain = 1e9;

    for (int level=0;; ++level) {
        std::cerr<<"--- Phase "<<level<<": n="<<g.n
                 <<", m="<<g.m<<" ---\n";

        bool doColor = (g.n > MIN_RECOLOR_VERTS
                        && lastPhaseGain >= TH_COLORED);

        std::vector<std::vector<int>> colorSets;
        if (doColor) {
            colorSets = createColorSets(g);
        } else {
            colorSets.clear();
            colorSets.push_back(std::vector<int>(g.n));
            std::iota(colorSets[0].begin(), colorSets[0].end(), 0);
        }

        // init communities
        std::vector<int> C(g.n);
        std::iota(C.begin(), C.end(), 0);

        // pick threshold
        double th = doColor ? TH_COLORED : TH_FINE;

        bool moved=true;
        int iter=0;
        while (moved && iter<100) {
            moved = localMovePhaseWithIsolateSets(g, C, colorSets, th);
            ++iter;
        }

        // compute this phase’s gain
        double Qpost = computeModularity(g,C);
        double Qprev = (H.partitions.empty()
                        ? computeModularity(g,C)  // same
                        : computeModularity(g,C));
        lastPhaseGain = Qpost - Qprev;

        // record partition on original vertices
        std::vector<int> part(input_g.n);
        for (int i=0; i<input_g.n; ++i)
            part[i] = C[hmap[i]];
        H.partitions.push_back(part);

        // aggregate
        std::unordered_map<int,int> c2n;
        Graph next = aggregateGraphFast(g, C, c2n);
        if (next.n == g.n) break;

        // update map
        for (int i=0; i<input_g.n; ++i)
            hmap[i] = c2n[hmap[i]];
        g = std::move(next);
    }

    if (H.partitions.empty()) {
        std::vector<int> triv(input_g.n);
        std::iota(triv.begin(), triv.end(), 0);
        H.partitions.push_back(triv);
    }
}
