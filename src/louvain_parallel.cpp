#include "louvain_parallel.h"

#include <algorithm>
#include <random>
#include <numeric>
#include <iostream>
#include <atomic>
#include <cmath>
#include <unordered_set>

// Partition the graph into roughly equal parts for parallel processing
std::vector<GraphPartition> partitionGraph(const Graph &g, int numPartitions) {
    int n = g.n;
    std::vector<GraphPartition> partitions(numPartitions);
    
    // Simple partitioning strategy: divide nodes approximately equally
    int nodesPerPartition = n / numPartitions;
    int remainingNodes = n % numPartitions;
    
    int currentNode = 0;
    for (int p = 0; p < numPartitions; p++) {
        // Calculate how many nodes this partition should get
        int partitionSize = nodesPerPartition + (p < remainingNodes ? 1 : 0);
        
        // Assign nodes to this partition
        for (int i = 0; i < partitionSize; i++) {
            if (currentNode < n) {
                partitions[p].nodeIds.push_back(currentNode);
                
                // Map global to local and local to global indices
                int localIdx = partitions[p].globalToLocal.size();
                partitions[p].globalToLocal[currentNode] = localIdx;
                partitions[p].localToGlobal[localIdx] = currentNode;
                
                currentNode++;
            }
        }
        
        // Initialize local adjacency list and degree vector
        int localSize = partitions[p].nodeIds.size();
        partitions[p].localAdj.resize(localSize);
        partitions[p].localDegree.resize(localSize, 0.0);
        
        // Populate local adjacency list and degree vector
        for (int i = 0; i < localSize; i++) {
            int globalNodeId = partitions[p].nodeIds[i];
            
            // Copy the degree
            partitions[p].localDegree[i] = g.degree[globalNodeId];
            
            // Copy adjacency list
            for (const auto &edge : g.adj[globalNodeId]) {
                partitions[p].localAdj[i].push_back(edge);
            }
        }
    }
    
    return partitions;
}

// Helper function to perform one pass of local moves on a partition
static bool parallelLocalMovePass(const Graph &g, const GraphPartition &partition,
                                 std::vector<int> &community, std::vector<double> &comm_tot,
                                 std::vector<int> &localMoves) {
    bool moved = false;
    double m_total = 2.0 * g.m;
    
    // Process nodes in this partition
    for (size_t idx = 0; idx < partition.nodeIds.size(); idx++) {
        int i = partition.nodeIds[idx];
        int current_comm = community[i];
        double k_i = g.degree[i];
        
        // Compute weight sum from node i to each neighbor community
        std::unordered_map<int, double> neighCommWeight;
        for (const auto &edge : g.adj[i]) {
            int nb_comm = community[edge.neighbor];
            neighCommWeight[nb_comm] += edge.weight;
        }
        
        // Thread-safe update of comm_tot using atomic operations
        #pragma omp atomic update
        comm_tot[current_comm] -= k_i;
        
        double best_gain = 0.0;
        int best_comm = current_comm;
        
        // Evaluate gain for each neighboring community
        for (const auto &entry : neighCommWeight) {
            int c = entry.first;
            double k_i_in = entry.second;
            
            // Read the current community total (safe for reading)
            double current_comm_tot;
            #pragma omp atomic read
            current_comm_tot = comm_tot[c];
            
            double gain = k_i_in - (current_comm_tot * k_i) / m_total;
            if (gain > best_gain) {
                best_gain = gain;
                best_comm = c;
            }
        }
        
        if (best_comm == -1)
            best_comm = current_comm;
        
        // Store local move decision
        localMoves[idx] = best_comm;
        
        if (best_comm != current_comm)
            moved = true;
    }
    
    return moved;
}

// Perform local move phase on a partition
std::vector<int> parallelLocalMovePhase(const Graph &g, const GraphPartition &partition, 
                                       const std::vector<int> &initialCommunities) {
    std::vector<int> community = initialCommunities;
    std::vector<double> comm_tot(g.n, 0.0);
    
    // Initialize community totals
    for (int i = 0; i < g.n; ++i) {
        comm_tot[community[i]] += g.degree[i];
    }
    
    double currentQ = computeModularity(g, community);
    bool improvement = true;
    int iter = 0;
    
    // Vector to store local move decisions
    std::vector<int> localMoves(partition.nodeIds.size());
    
    while (improvement && iter < 10) {  // Limit iterations for safety
        iter++;
        
        // Calculate local moves without applying them yet
        bool moved = parallelLocalMovePass(g, partition, community, comm_tot, localMoves);
        
        // Apply the local moves
        #pragma omp critical
        {
            for (size_t idx = 0; idx < partition.nodeIds.size(); idx++) {
                int i = partition.nodeIds[idx];
                int new_comm = localMoves[idx];
                
                if (community[i] != new_comm) {
                    community[i] = new_comm;
                    comm_tot[new_comm] += g.degree[i];
                }
            }
        }
        
        double newQ = computeModularity(g, community);
        
        #pragma omp critical
        {
            std::cout << "Thread " << omp_get_thread_num() 
                      << " - Local move pass " << iter 
                      << " on partition: modularity = " << newQ << "\n";
        }
        
        if (std::fabs(newQ - currentQ) < 1e-6)
            break;
        
        currentQ = newQ;
        improvement = moved;
    }
    
    return community;
}

void louvainParallel(const Graph &g, Hierarchy &hierarchy, int numThreads) {
    // Set number of OpenMP threads
    omp_set_num_threads(numThreads);
    std::cout << "Running parallel Louvain with " << numThreads << " threads\n";
    
    // Create a copy of the original graph for recursive calls
    Graph current_g = g;
    
    // Initialize full hierarchy mapping
    std::vector<int> hierarchy_mapping(g.n);
    for (int i = 0; i < g.n; i++) {
        hierarchy_mapping[i] = i;
    }
    
    int level = 0;
    bool improvement = true;
    
    while (improvement && level < 10) {  // Limit levels for safety
        // Initialize community assignments (each node in its own community)
        std::vector<int> current_communities(current_g.n);
        for (int i = 0; i < current_g.n; i++) {
            current_communities[i] = i;
        }
        
        double currentQ = computeModularity(current_g, current_communities);
        double prevQ = currentQ - 1.0; // Ensure first iteration runs
        int iter = 0;
        
        // Perform local move phase with iterative refinement
        while (currentQ > prevQ + 1e-6 && iter < 20) {  // Limit iterations and check for improvement
            iter++;
            prevQ = currentQ;
            
            // Use a two-phase approach: compute moves first, then apply them
            std::vector<int> best_communities(current_g.n, -1);
            
            // Compute community weights (totals) - shared by all threads
            std::vector<double> comm_tot(current_g.n, 0.0);
            for (int i = 0; i < current_g.n; i++) {
                comm_tot[current_communities[i]] += current_g.degree[i];
            }
            
            // Phase 1: Compute best moves for all nodes in parallel
            #pragma omp parallel
            {
                #pragma omp for schedule(dynamic, 64)
                for (int i = 0; i < current_g.n; i++) {
                    int current_comm = current_communities[i];
                    double k_i = current_g.degree[i];
                    double m_total = 2.0 * current_g.m;
                    
                    // Skip isolated nodes (no neighbors)
                    if (current_g.adj[i].empty()) {
                        best_communities[i] = current_comm;
                        continue;
                    }
                    
                    // Compute weight sum from node i to each neighboring community
                    std::unordered_map<int, double> neighCommWeight;
                    for (const auto &edge : current_g.adj[i]) {
                        int nb_comm = current_communities[edge.neighbor];
                        neighCommWeight[nb_comm] += edge.weight;
                    }
                    
                    // Find best community
                    double best_gain = 0.0;
                    int best_comm = current_comm;
                    
                    // Remove contribution from current community
                    double self_loop_weight = 0.0;
                    for (const auto &edge : current_g.adj[i]) {
                        if (edge.neighbor == i) {
                            self_loop_weight += edge.weight;
                        }
                    }
                    
                    // Calculate gain for each neighboring community
                    for (const auto &entry : neighCommWeight) {
                        int c = entry.first;
                        double k_i_in = entry.second;
                        
                        // Skip current community in calculation
                        if (c == current_comm) continue;
                        
                        // Thread-safe read of community total
                        double c_tot;
                        #pragma omp atomic read
                        c_tot = comm_tot[c];
                        
                        // Calculate modularity gain
                        double comm_tot_without_i = (c == current_comm) ? 
                            (c_tot - k_i) : c_tot;
                        
                        double gain = k_i_in - (comm_tot_without_i * k_i) / m_total;
                        
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_comm = c;
                        }
                    }
                    
                    // Store best community
                    best_communities[i] = best_comm;
                }
            }
            
            // Phase 2: Apply moves sequentially to avoid race conditions
            bool any_moved = false;
            for (int i = 0; i < current_g.n; i++) {
                int new_comm = best_communities[i];
                int old_comm = current_communities[i];
                
                if (new_comm != old_comm && new_comm != -1) {
                    current_communities[i] = new_comm;
                    any_moved = true;
                }
            }
            
            // Calculate new modularity
            currentQ = computeModularity(current_g, current_communities);
            std::cout << "Level " << level << " - Pass " << iter 
                      << ": modularity = " << currentQ << "\n";
            
            // Stop if no nodes moved
            if (!any_moved) break;
        }
        
        // Count unique communities
        std::unordered_set<int> uniqueComm(current_communities.begin(), current_communities.end());
        std::cout << "Level " << level << " final modularity: " << currentQ << ", "
                  << "communities = " << uniqueComm.size() << "\n";
        
        // If we've reached a point where each node is alone, stop
        if (uniqueComm.size() == (size_t)current_g.n || uniqueComm.size() <= 1) {
            improvement = false;
            break;
        }
        
        // Store current level mapping
        std::vector<int> level_mapping(g.n);
        for (int i = 0; i < g.n; i++) {
            if (i < hierarchy_mapping.size()) {
                int old_comm = hierarchy_mapping[i];
                if (old_comm < (int)current_communities.size()) {
                    level_mapping[i] = current_communities[old_comm];
                } else {
                    level_mapping[i] = i;
                }
            } else {
                level_mapping[i] = i;
            }
        }
        
        // Add this level's partition to the hierarchy
        hierarchy.partitions.push_back(level_mapping);
        
        // Aggregate graph
        std::unordered_map<int, int> commMap;
        Graph g_new = aggregateGraph(current_g, current_communities, commMap);
        
        // Update hierarchy mapping
        std::vector<int> new_hierarchy_mapping(g.n);
        for (int i = 0; i < g.n; i++) {
            if (i < hierarchy_mapping.size()) {
                int old_comm = hierarchy_mapping[i];
                if (old_comm < (int)current_communities.size() && 
                    current_communities[old_comm] < (int)commMap.size()) {
                    new_hierarchy_mapping[i] = commMap[current_communities[old_comm]];
                } else {
                    new_hierarchy_mapping[i] = i;
                }
            } else {
                new_hierarchy_mapping[i] = i;
            }
        }
        hierarchy_mapping = new_hierarchy_mapping;
        
        // If graph couldn't be further aggregated, stop
        if (g_new.n >= current_g.n || g_new.n <= 1) {
            improvement = false;
            break;
        }
        
        // Update current graph for next iteration
        current_g = g_new;
        level++;
    }
    
    // If no partitions were added, add trivial partition (each node in its own community)
    if (hierarchy.partitions.empty()) {
        std::vector<int> trivial_partition(g.n);
        for (int i = 0; i < g.n; i++) {
            trivial_partition[i] = i;
        }
        hierarchy.partitions.push_back(trivial_partition);
    }
}