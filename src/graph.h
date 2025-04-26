#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>

// A structure for a (weighted) edge.
struct Edge {
    int neighbor;    // Index of the neighbor node.
    double weight;   // Weight of the edge.
};

// Graph structure: 
//  - n = number of nodes
//  - m = total weight of undirected edges
//  - adj = adjacency list
//  - degree[i] = weighted degree of node i
struct Graph {
    int n;
    double m;
    std::vector<std::vector<Edge>> adj;
    std::vector<double> degree;
};

// Load graph from a text file of the form:
//   <n> <m>
//   u v [weight]
// (skips comment lines starting with # or %)
Graph readGraph(const std::string &filename);

// Compute modularity of a partition:
//  partition[i] = community of node i
double computeModularity(const Graph &g, const std::vector<int> &partition);


#endif // GRAPH_H
