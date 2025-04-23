#include "louvain_parallel_vfc.h"

#include <algorithm>
#include <random>
#include <numeric>
#include <iostream>
#include <atomic>
#include <cmath>
#include <unordered_set>
#include <map>
#include <vector>
#include <queue>
#include <chrono>
#include <mutex>

// Create isolate sets for conflict-free parallel processing
// Isolate sets contain vertices that don't share neighbors
std::vector<std::vector<int>> createIsolateSets(const Graph &g) {
    // Simplified version that uses distance-1 coloring
    std::vector<int> colors(g.n, -1); // -1 means uncolored
    std::vector<std::vector<int>> isolateSets;
    
    // For deterministic results
    std::vector<int> vertices(g.n);
    std::iota(vertices.begin(), vertices.end(), 0);
    
    // Simple greedy coloring algorithm
    for (int v : vertices) {
        // Find colors of neighbors
        std::vector<bool> used_colors(g.n, false);
        for (const auto& edge : g.adj[v]) {
            int neighbor = edge.neighbor;
            if (colors[neighbor] >= 0) {
                used_colors[colors[neighbor]] = true;
            }
        }
        
        // Find the smallest unused color
        int color = 0;
        while (color < g.n && used_colors[color]) {
            color++;
        }
        
        colors[v] = color;
    }
    
    // Group vertices by color
    int max_color = *std::max_element(colors.begin(), colors.end());
    isolateSets.resize(max_color + 1);
    
    for (int i = 0; i < g.n; i++) {
        isolateSets[colors[i]].push_back(i);
    }
    
    // Remove empty color sets
    isolateSets.erase(
        std::remove_if(isolateSets.begin(), isolateSets.end(), 
                      [](const std::vector<int>& set) { return set.empty(); }),
        isolateSets.end()
    );
    
    std::cout << "Created " << isolateSets.size() << " color sets\n";
    return isolateSets;
}

// Optimized Vertex Following implementation
Graph applyVertexFollowing(const Graph &g, std::vector<int> &mapping) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    mapping.resize(g.n);
    std::iota(mapping.begin(), mapping.end(), 0);
    
    // Identify single-degree vertices
    std::vector<bool> is_single_degree(g.n, false);
    std::vector<int> neighbor_of_single(g.n, -1);
    int single_degree_count = 0;
    
    #pragma omp parallel for reduction(+:single_degree_count)
    for (int i = 0; i < g.n; i++) {
        if (g.adj[i].size() == 1) {
            is_single_degree[i] = true;
            neighbor_of_single[i] = g.adj[i][0].neighbor;
            single_degree_count++;
        }
    }
    
    // If there are no single degree vertices, just return the original graph
    if (single_degree_count == 0) {
        return g;
    }
    
    // Create a mapping from original vertices to new vertices
    int new_vertex_count = 0;
    std::vector<int> old_to_new(g.n, -1);
    
    for (int i = 0; i < g.n; i++) {
        if (is_single_degree[i]) {
            // Map this vertex to its neighbor
            mapping[i] = neighbor_of_single[i];
        } else {
            // This vertex remains
            old_to_new[i] = new_vertex_count++;
        }
    }
    
    // Update the mapping to use new vertex indices
    for (int i = 0; i < g.n; i++) {
        if (is_single_degree[i]) {
            int target = mapping[i];
            // Follow the chain until we reach a non-single-degree vertex
            while (is_single_degree[target]) {
                target = mapping[target];
            }
            mapping[i] = old_to_new[target];
        } else {
            mapping[i] = old_to_new[i];
        }
    }
    
    // Create the new graph
    Graph new_g;
    new_g.n = new_vertex_count;
    new_g.adj.resize(new_g.n);
    new_g.degree.assign(new_g.n, 0.0);
    
    // Use thread-local storage for edge accumulation to reduce synchronization
    std::vector<std::vector<std::map<int, double>>> thread_local_edges;
    int num_threads = omp_get_max_threads();
    thread_local_edges.resize(num_threads);
    for (int t = 0; t < num_threads; t++) {
        thread_local_edges[t].resize(new_g.n);
    }
    
    // Add edges to the new graph in parallel
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < g.n; i++) {
            if (is_single_degree[i]) continue;  // Skip single-degree vertices
            
            int new_i = old_to_new[i];
            
            for (const auto &edge : g.adj[i]) {
                int j = edge.neighbor;
                double weight = edge.weight;
                
                int new_j = mapping[j];
                
                if (new_i != new_j) {
                    // Add to thread-local storage
                    thread_local_edges[thread_id][new_i][new_j] += weight;
                } else {
                    // Self-loop
                    thread_local_edges[thread_id][new_i][new_i] += weight;
                }
            }
        }
    }
    
    // Merge thread-local edges
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < new_g.n; i++) {
        // Merge edges from all threads
        std::map<int, double> merged_edges;
        for (int t = 0; t < num_threads; t++) {
            for (const auto &entry : thread_local_edges[t][i]) {
                merged_edges[entry.first] += entry.second;
            }
        }
        
        // Create final edge list and calculate degrees
        for (const auto &entry : merged_edges) {
            int j = entry.first;
            double weight = entry.second;
            
            new_g.adj[i].push_back({j, weight});
            new_g.degree[i] += weight;
            
            // Count edge once for total weight (if it's not a self-loop)
            if (i <= j && i != j) {
                #pragma omp atomic
                new_g.m += weight / 2.0;
            } else if (i == j) {
                // Handle self-loops specially
                #pragma omp atomic
                new_g.m += weight / 2.0;
            }
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    
    std::cout << "Vertex Following reduced graph from " << g.n << " to " << new_g.n 
              << " vertices (" << (100.0 * (g.n - new_g.n) / g.n) << "% reduction) in "
              << duration << " ms\n";
    
    return new_g;
}

// Perform a local move phase using isolate sets for conflict-free parallelism
bool localMovePhaseWithIsolateSets(Graph &g, std::vector<int> &communities, 
                                 const std::vector<std::vector<int>> &isolateSets,
                                 double threshold) {
    double m_total = 2.0 * g.m;
    bool any_moved = false;
    double initial_modularity = computeModularity(g, communities);
    std::cout << "Starting local move phase, initial modularity: " << initial_modularity << std::endl;
    
    // Process each isolate set
    for (size_t set_idx = 0; set_idx < isolateSets.size(); set_idx++) {
        const auto &isolateSet = isolateSets[set_idx];
        
        // Recalculate community totals for this iteration
        std::vector<double> comm_tot(g.n, 0.0);
        for (int i = 0; i < g.n; i++) {
            comm_tot[communities[i]] += g.degree[i];
        }
        
        // Process vertices in this isolate set in parallel
        bool local_any_moved = false;
        
        #pragma omp parallel
        {
            bool thread_moved = false;
            
            #pragma omp for schedule(dynamic, 64)
            for (size_t idx = 0; idx < isolateSet.size(); idx++) {
                int i = isolateSet[idx];
                int current_comm = communities[i];
                double k_i = g.degree[i];
                
                // Skip isolated vertices (no neighbors)
                if (g.adj[i].empty()) continue;
                
                // Compute neighboring communities and their weights
                std::unordered_map<int, double> neighCommWeight;
                for (const auto &edge : g.adj[i]) {
                    int nb_comm = communities[edge.neighbor];
                    neighCommWeight[nb_comm] += edge.weight;
                }
                
                // Temporarily remove node from its community
                double self_loop_weight = 0.0;
                for (const auto &edge : g.adj[i]) {
                    if (edge.neighbor == i) {
                        self_loop_weight += edge.weight;
                        break;
                    }
                }
                
                // Find best community
                double best_gain = 0.0;
                int best_comm = current_comm;
                
                // Evaluate gain for each neighboring community
                for (const auto &entry : neighCommWeight) {
                    int c = entry.first;
                    double weight_to_c = entry.second;
                    
                    // Skip current community in gain calculation
                    if (c == current_comm) continue;
                    
                    // Calculate modularity gain
                    double gain = weight_to_c - (k_i * comm_tot[c]) / m_total;
                    
                    if (gain > best_gain || (gain == best_gain && c < best_comm)) {
                        best_gain = gain;
                        best_comm = c;
                    }
                }
                
                // Update community assignment if beneficial
                if (best_comm != current_comm) {
                    communities[i] = best_comm;
                    thread_moved = true;
                }
            }
            
            // Combine thread results
            if (thread_moved) {
                #pragma omp critical
                {
                    local_any_moved = true;
                }
            }
        }
        
        if (local_any_moved) {
            any_moved = true;
        }
        
        // Print progress for every few sets processed
        if (set_idx % 10 == 0 || set_idx == isolateSets.size() - 1) {
            std::cout << "  Processed " << (set_idx + 1) << "/" << isolateSets.size() 
                      << " isolate sets\r" << std::flush;
        }
    }
    std::cout << std::endl;
    
    double final_modularity = computeModularity(g, communities);
    double gain = final_modularity - initial_modularity;
    
    std::cout << "Local move phase completed, final modularity: " << final_modularity 
              << ", gain: " << gain << std::endl;
    
    // Return whether the phase should continue
    return any_moved && (gain > threshold * std::abs(initial_modularity));
}

// Create an aggregated graph based on communities
Graph aggregateGraphFast(const Graph &g, const std::vector<int> &communities, 
                       std::unordered_map<int, int> &commMap) {
    commMap.clear();
    
    // Identify unique communities and map them to consecutive indices
    std::vector<int> unique_comms;
    for (int i = 0; i < g.n; i++) {
        if (commMap.find(communities[i]) == commMap.end()) {
            commMap[communities[i]] = commMap.size();
            unique_comms.push_back(communities[i]);
        }
    }
    
    // Create new graph
    Graph new_g;
    new_g.n = unique_comms.size();
    new_g.adj.resize(new_g.n);
    new_g.degree.assign(new_g.n, 0.0);
    new_g.m = 0.0;
    
    // Use thread-local storage to avoid synchronization
    int num_threads = omp_get_max_threads();
    std::vector<std::vector<std::map<int, double>>> thread_local_edges(num_threads);
    for (int t = 0; t < num_threads; t++) {
        thread_local_edges[t].resize(new_g.n);
    }
    
    // Process edges in parallel to create community graph
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < g.n; i++) {
            int comm_i = commMap[communities[i]];
            
            for (const auto &edge : g.adj[i]) {
                int j = edge.neighbor;
                double weight = edge.weight;
                int comm_j = commMap[communities[j]];
                
                // Add to thread-local storage
                thread_local_edges[thread_id][comm_i][comm_j] += weight;
            }
        }
    }
    
    // Merge thread-local results and build the final graph
    for (int i = 0; i < new_g.n; i++) {
        // Merge edges from all threads
        std::map<int, double> merged_edges;
        for (int t = 0; t < num_threads; t++) {
            for (const auto &entry : thread_local_edges[t][i]) {
                merged_edges[entry.first] += entry.second;
            }
        }
        
        // Create final edge list and calculate degrees
        for (const auto &entry : merged_edges) {
            int j = entry.first;
            double weight = entry.second;
            
            new_g.adj[i].push_back({j, weight});
            new_g.degree[i] += weight;
            
            // Count edge once for total weight
            if (i <= j) {
                new_g.m += weight / 2.0;
            }
        }
    }
    
    return new_g;
}

// Run the Louvain algorithm with optimized VF and isolate sets
void louvainParallelVFC(const Graph &input_g, Hierarchy &hierarchy, int numThreads) {
    // Set number of OpenMP threads
    omp_set_num_threads(numThreads);
    std::cout << "Running parallel Louvain with VF+Color using " 
              << numThreads << " threads\n";
    
    // Apply Vertex Following heuristic
    std::vector<int> vf_mapping;
    Graph g = input_g;
    
    // Skip vertex following for now to debug the main algorithm
    /*
    auto vf_start = std::chrono::high_resolution_clock::now();
    g = applyVertexFollowing(input_g, vf_mapping);
    auto vf_end = std::chrono::high_resolution_clock::now();
    auto vf_duration = std::chrono::duration_cast<std::chrono::milliseconds>(vf_end - vf_start).count();
    */
    
    // Create a map to track the full hierarchy from original vertices
    std::vector<int> hierarchy_mapping(input_g.n);
    for (int i = 0; i < input_g.n; i++) {
        hierarchy_mapping[i] = i; //vf_mapping[i];
    }
    
    // Run the algorithm
    Graph current_g = g;
    int level = 0;
    bool improvement = true;
    
    // Modularity gain thresholds
    const double THRESHOLD = 1e-4;  // Single threshold for simplicity during debugging
    
    while (improvement && level < 3) {  // Limit levels for debugging
        std::cout << "\n--- Starting Phase " << level << " ---\n";
        
        // Create color sets for parallel processing
        auto color_start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<int>> isolateSets = createIsolateSets(current_g);
        auto color_end = std::chrono::high_resolution_clock::now();
        auto color_duration = std::chrono::duration_cast<std::chrono::milliseconds>(color_end - color_start).count();
        
        std::cout << "Created " << isolateSets.size() << " color sets in " 
                  << color_duration << " ms\n";
        
        // Initialize community assignments (each node in its own community)
        std::vector<int> current_communities(current_g.n);
        for (int i = 0; i < current_g.n; i++) {
            current_communities[i] = i;
        }
        
        double currentQ = computeModularity(current_g, current_communities);
        double prevQ = currentQ - 1.0;  // Ensure first iteration runs
        int iter = 0;
        bool phase_improvement = true;
        
        std::cout << "Initial modularity: " << currentQ << "\n";
        
        // Perform iterations until convergence
        while (phase_improvement && iter < 10) {  // Limit iterations for debugging
            iter++;
            prevQ = currentQ;
            
            // Perform local move phase using isolate sets
            phase_improvement = localMovePhaseWithIsolateSets(
                current_g, current_communities, isolateSets, THRESHOLD);
            
            // Calculate new modularity
            currentQ = computeModularity(current_g, current_communities);
            std::cout << "Phase " << level << " - Iteration " << iter 
                      << ": modularity = " << currentQ
                      << ", gain = " << (currentQ - prevQ) << "\n";
            
            // Check for convergence
            if (std::fabs(currentQ - prevQ) < THRESHOLD * std::abs(prevQ) || !phase_improvement) {
                break;
            }
        }
        
        // Count unique communities
        std::unordered_set<int> uniqueComm(current_communities.begin(), current_communities.end());
        std::cout << "Phase " << level << " final: modularity = " << currentQ 
                  << ", communities = " << uniqueComm.size() << " (from " << current_g.n << " vertices)\n";
        
        // If minimal or no community reduction, stop
        if (uniqueComm.size() >= current_g.n - 1 || uniqueComm.size() <= 1) {
            std::cout << "No significant reduction in communities, stopping.\n";
            improvement = false;
            break;
        }
        
        // Update the hierarchy mapping
        std::vector<int> level_mapping(input_g.n);
        for (int i = 0; i < input_g.n; i++) {
            int old_comm = hierarchy_mapping[i];
            if (old_comm < (int)current_communities.size()) {
                level_mapping[i] = current_communities[old_comm];
            } else {
                level_mapping[i] = old_comm;
            }
        }
        
        // Add this level's partition to the hierarchy
        hierarchy.partitions.push_back(level_mapping);
        hierarchy_mapping = level_mapping;
        
        // Prepare for aggregation
        std::unordered_map<int, int> commMap;
        for (int i = 0; i < current_g.n; i++) {
            int comm = current_communities[i];
            if (commMap.find(comm) == commMap.end()) {
                commMap[comm] = commMap.size();
            }
        }
        
        // Create new aggregated graph
        Graph next_g;
        next_g.n = commMap.size();
        next_g.adj.resize(next_g.n);
        next_g.degree.assign(next_g.n, 0.0);
        next_g.m = 0.0;
        
        // Build the aggregated graph
        for (int i = 0; i < current_g.n; i++) {
            int comm_i = commMap[current_communities[i]];
            
            for (const auto &edge : current_g.adj[i]) {
                int j = edge.neighbor;
                double weight = edge.weight;
                int comm_j = commMap[current_communities[j]];
                
                // Add edge to the new graph
                bool found = false;
                for (auto &e : next_g.adj[comm_i]) {
                    if (e.neighbor == comm_j) {
                        e.weight += weight;
                        found = true;
                        break;
                    }
                }
                
                if (!found) {
                    next_g.adj[comm_i].push_back({comm_j, weight});
                }
                
                // Update degree
                next_g.degree[comm_i] += weight;
                
                // Count each edge once for total weight
                if (i <= j) {
                    next_g.m += weight / 2.0;
                }
            }
        }
        
        std::cout << "Aggregated graph: " << current_g.n << " -> " << next_g.n << " vertices\n";
        
        // Update current graph for next phase
        current_g = next_g;
        level++;
    }
    
    // If no partitions were added, add trivial partition
    if (hierarchy.partitions.empty()) {
        std::vector<int> trivial_partition(input_g.n);
        for (int i = 0; i < input_g.n; i++) {
            trivial_partition[i] = i;
        }
        hierarchy.partitions.push_back(trivial_partition);
    }
}