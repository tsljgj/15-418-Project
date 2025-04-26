#include "graph.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <numeric>

Graph readGraph(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) { std::cerr<<"Cannot open "<<filename<<"\n"; exit(1); }
    std::string line;
    while (std::getline(in,line)) {
        if (line.empty()||line[0]=='#'||line[0]=='%') continue;
        break;
    }
    int n, m_dummy;
    { std::istringstream iss(line); iss >> n >> m_dummy; }
    Graph g;
    g.n = n;
    g.adj.assign(n,{});
    g.degree.assign(n,0.0);

    int u,v; double w;
    while (std::getline(in,line)) {
        if (line.empty()||line[0]=='#'||line[0]=='%') continue;
        std::istringstream iss(line);
        if (!(iss>>u>>v)) continue;
        if (!(iss>>w)) w = 1.0;
        if (u<0||u>=n||v<0||v>=n||u==v) continue;
        g.adj[u].push_back({v,w});
        g.adj[v].push_back({u,w});
        g.degree[u] += w;
        g.degree[v] += w;
    }
    double tot = std::accumulate(g.degree.begin(), g.degree.end(), 0.0);
    g.m = tot/2.0;
    return g;
}

double computeModularity(const Graph &g, const std::vector<int> &comm) {
    int n = g.n;
    double two_m = 2.0 * g.m;
    std::unordered_map<int,double> tot, in;
    tot.reserve(n); in.reserve(n);
    for (int i = 0; i < n; i++) {
        tot[comm[i]] += g.degree[i];
        for (auto &e : g.adj[i]) {
            if (comm[i]==comm[e.neighbor]) in[comm[i]] += e.weight;
        }
    }
    double Q = 0;
    for (auto &p : tot) {
        int c = p.first;
        double t = p.second;
        double iw = in[c]/2.0;
        Q += (iw/g.m) - (t*t)/(two_m*two_m);
    }
    return Q;
}

