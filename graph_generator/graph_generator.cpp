#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_set>
#include <string>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <filesystem>  // new include

// Struct to represent an edge
struct Edge {
    int source;
    int target;
    double weight;
};

// Function to write the graph to a file in the format expected by the Louvain algorithm
void writeGraphToFile(const std::string& filename, int numNodes, const std::vector<Edge>& edges) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    // Write header: number of nodes and edges
    file << numNodes << " " << edges.size() << std::endl;
    
    // Write edges
    for (const auto& edge : edges) {
        file << edge.source << " " << edge.target << " " << edge.weight << std::endl;
    }
    
    file.close();
    std::cout << "Graph written to " << filename << std::endl;
}

// Generate a random graph with n nodes and approximately m edges
std::vector<Edge> generateRandomGraph(int n, int m, double minWeight, double maxWeight, 
                                     std::mt19937& gen) {
    std::vector<Edge> edges;
    std::uniform_real_distribution<double> weightDist(minWeight, maxWeight);
    std::uniform_int_distribution<int> nodeDist(0, n-1);
    
    // Using a set to avoid duplicate edges
    std::unordered_set<std::string> edgeSet;
    
    while (edges.size() < static_cast<size_t>(m)) {
        int source = nodeDist(gen);
        int target = nodeDist(gen);
        
        // Avoid self-loops
        if (source == target) continue;
        
        // Make sure we have source < target to avoid duplicates (undirected graph)
        if (source > target) std::swap(source, target);
        
        // Create a unique string representation of the edge
        std::string edgeStr = std::to_string(source) + "-" + std::to_string(target);
        
        // Add edge if it doesn't already exist
        if (edgeSet.find(edgeStr) == edgeSet.end()) {
            double weight = weightDist(gen);
            edges.push_back({source, target, weight});
            edgeSet.insert(edgeStr);
        }
    }
    
    return edges;
}

// Generate a large community graph efficiently - optimized for scale
std::vector<Edge> generateCommunityGraph(int numNodes, int numCommunities, double pIntra, double pInter,
                                           double minWeight, double maxWeight, std::mt19937& gen) {
    std::vector<Edge> edges;
    std::uniform_real_distribution<double> weightDist(minWeight, maxWeight);
    
    // Assign nodes to communities
    std::vector<int> communities(numNodes);
    int nodesPerCommunity = numNodes / numCommunities;
    for (int i = 0; i < numNodes; i++) {
        communities[i] = std::min(i / nodesPerCommunity, numCommunities - 1);
    }
    
    // Calculate community sizes
    std::vector<int> communitySizes(numCommunities, 0);
    for (int i = 0; i < numNodes; i++) {
        communitySizes[communities[i]]++;
    }
    
    // For each community, generate intra-community edges
    for (int c = 0; c < numCommunities; c++) {
        int communitySize = communitySizes[c];
        
        // For very large communities, use expected number of edges instead of checking all pairs
        if (communitySize > 10000) {
            // Expected number of edges in this community
            long long expectedIntraEdges = (long long)((communitySize * (communitySize - 1) / 2) * pIntra);
            std::cout << "Community " << c << " - Size: " << communitySize 
                      << " - Expected intra edges: " << expectedIntraEdges << std::endl;
            
            // Generate these edges randomly
            std::uniform_int_distribution<int> nodeDist(0, communitySize - 1);
            
            // Find the node range for this community
            int startNode = c * nodesPerCommunity;
            int endNode = (c == numCommunities - 1) ? numNodes - 1 : (c + 1) * nodesPerCommunity - 1;
            
            // Create a set to avoid duplicate edges
            std::unordered_set<std::string> edgeSet;
            
            // Progress tracking
            long long targetEdges = expectedIntraEdges;
            long long progressStep = targetEdges / 10; // Report progress every 10%
            if (progressStep == 0) progressStep = 1;
            
            std::cout << "Generating intra-community edges for community " << c << "..." << std::endl;
            for (long long e = 0; e < targetEdges && edgeSet.size() < targetEdges; e++) {
                // Progress reporting
                if (e % progressStep == 0) {
                    int percentComplete = (e * 100) / targetEdges;
                    std::cout << "  " << percentComplete << "% complete (" << e << "/" << targetEdges << " edges)" << std::endl;
                }
                
                // Generate random node pair within community
                int nodeOffset1 = nodeDist(gen);
                int nodeOffset2 = nodeDist(gen);
                
                if (nodeOffset1 == nodeOffset2) continue; // Skip self-loops
                
                int source = startNode + nodeOffset1;
                int target = startNode + nodeOffset2;
                
                if (source > target) std::swap(source, target);
                
                std::string edgeStr = std::to_string(source) + "-" + std::to_string(target);
                
                if (edgeSet.find(edgeStr) == edgeSet.end()) {
                    double weight = weightDist(gen);
                    edges.push_back({source, target, weight});
                    edgeSet.insert(edgeStr);
                }
            }
        } 
        else {
            // For smaller communities, use the original approach
            int startNode = c * nodesPerCommunity;
            int endNode = (c == numCommunities - 1) ? numNodes - 1 : (c + 1) * nodesPerCommunity - 1;
            
            std::uniform_real_distribution<double> randProb(0.0, 1.0);
            
            for (int i = startNode; i < endNode; i++) {
                for (int j = i + 1; j <= endNode; j++) {
                    if (randProb(gen) < pIntra) {
                        double weight = weightDist(gen);
                        edges.push_back({i, j, weight});
                    }
                }
            }
        }
    }
    
    // Generate inter-community edges - use sampling approach for efficiency
    if (numCommunities > 1) {
        std::cout << "Generating inter-community edges..." << std::endl;
        
        // For large graphs, we use a sampling approach to generate some inter-community edges
        long long totalPossibleInterEdges = 0;
        for (int c1 = 0; c1 < numCommunities; c1++) {
            for (int c2 = c1 + 1; c2 < numCommunities; c2++) {
                totalPossibleInterEdges += (long long)communitySizes[c1] * communitySizes[c2];
            }
        }
        
        // Either generate expected number of edges or use a cap for very large graphs
        long long targetInterEdges = std::min(
            (long long)(totalPossibleInterEdges * pInter),
            (long long)(numNodes * 10) // Cap to ensure we don't try to generate too many edges
        );
        
        std::cout << "Target inter-community edges: " << targetInterEdges << std::endl;
        
        // Create a set to track generated edges
        std::unordered_set<std::string> interEdgeSet;
        
        // Progress tracking
        long long progressStep = targetInterEdges / 10;
        if (progressStep == 0) progressStep = 1;
        
        // Uniform distributions for selecting communities and nodes
        std::uniform_int_distribution<int> commDist(0, numCommunities - 1);
        
        for (long long e = 0; e < targetInterEdges * 2 && interEdgeSet.size() < targetInterEdges; e++) {
            // Progress reporting
            if (e % progressStep == 0) {
                int percentComplete = (interEdgeSet.size() * 100) / targetInterEdges;
                std::cout << "  " << percentComplete << "% complete (" 
                          << interEdgeSet.size() << "/" << targetInterEdges << " edges)" << std::endl;
            }
            
            // Select two different communities
            int c1 = commDist(gen);
            int c2 = commDist(gen);
            if (c1 == c2) continue;
            
            // Get community ranges
            int c1Start = c1 * nodesPerCommunity;
            int c1End = (c1 == numCommunities - 1) ? numNodes - 1 : (c1 + 1) * nodesPerCommunity - 1;
            
            int c2Start = c2 * nodesPerCommunity;
            int c2End = (c2 == numCommunities - 1) ? numNodes - 1 : (c2 + 1) * nodesPerCommunity - 1;
            
            // Select random nodes from each community
            std::uniform_int_distribution<int> node1Dist(c1Start, c1End);
            std::uniform_int_distribution<int> node2Dist(c2Start, c2End);
            
            int source = node1Dist(gen);
            int target = node2Dist(gen);
            
            if (source > target) std::swap(source, target);
            
            std::string edgeStr = std::to_string(source) + "-" + std::to_string(target);
            
            if (interEdgeSet.find(edgeStr) == interEdgeSet.end()) {
                double weight = weightDist(gen);
                edges.push_back({source, target, weight});
                interEdgeSet.insert(edgeStr);
            }
        }
    }
    
    // Ensure the graph is connected by adding edges between consecutive communities if needed
    for (int c = 0; c < numCommunities - 1; c++) {
        int sourceNode = c * nodesPerCommunity;
        int targetNode = (c + 1) * nodesPerCommunity;
        
        // Add a connecting edge with random weight
        double weight = weightDist(gen);
        edges.push_back({sourceNode, targetNode, weight});
    }
    
    return edges;
}

// Generate a preferential attachment graph (Barabási–Albert model)
std::vector<Edge> generatePreferentialAttachmentGraph(int numNodes, int m0, int m, 
                                                     double minWeight, double maxWeight,
                                                     std::mt19937& gen) {
    if (numNodes <= m0) {
        throw std::invalid_argument("Number of nodes must be greater than initial nodes (m0)");
    }
    if (m > m0) {
        throw std::invalid_argument("Edges per new node (m) must be less than or equal to initial nodes (m0)");
    }
    
    std::vector<Edge> edges;
    std::uniform_real_distribution<double> weightDist(minWeight, maxWeight);
    
    // Create a complete graph for the initial m0 nodes
    for (int i = 0; i < m0; i++) {
        for (int j = i + 1; j < m0; j++) {
            double weight = weightDist(gen);
            edges.push_back({i, j, weight});
        }
    }
    
    // Vector to track the degree of each node (for preferential attachment)
    std::vector<int> degrees(numNodes, 0);
    for (const auto& edge : edges) {
        degrees[edge.source]++;
        degrees[edge.target]++;
    }
    
    // Add remaining nodes with preferential attachment
    for (int i = m0; i < numNodes; i++) {
        // Calculate sum of all degrees for probability calculation
        int sumDegrees = 0;
        for (int j = 0; j < i; j++) {
            sumDegrees += degrees[j];
        }
        
        // Add m edges from the new node to existing nodes
        std::unordered_set<int> connectedNodes;
        while (connectedNodes.size() < static_cast<size_t>(m)) {
            // Choose an existing node with probability proportional to its degree
            std::uniform_int_distribution<int> dist(0, sumDegrees - 1);
            int r = dist(gen);
            int sum = 0;
            int target = -1;
            
            for (int j = 0; j < i; j++) {
                sum += degrees[j];
                if (r < sum) {
                    target = j;
                    break;
                }
            }
            
            if (target != -1 && connectedNodes.find(target) == connectedNodes.end()) {
                connectedNodes.insert(target);
                double weight = weightDist(gen);
                edges.push_back({target, i, weight});
                degrees[target]++;
                degrees[i]++;
            }
        }
    }
    
    return edges;
}

// Generate a small-world graph (Watts-Strogatz model)
std::vector<Edge> generateSmallWorldGraph(int numNodes, int k, double beta,
                                         double minWeight, double maxWeight,
                                         std::mt19937& gen) {
    if (k % 2 != 0) {
        throw std::invalid_argument("k must be even");
    }
    if (k >= numNodes) {
        throw std::invalid_argument("k must be less than numNodes");
    }
    
    std::vector<Edge> edges;
    std::uniform_real_distribution<double> weightDist(minWeight, maxWeight);
    std::uniform_real_distribution<double> randProb(0.0, 1.0);
    
    // Create a ring lattice
    for (int i = 0; i < numNodes; i++) {
        for (int j = 1; j <= k/2; j++) {
            int target = (i + j) % numNodes;
            double weight = weightDist(gen);
            edges.push_back({i, target, weight});
        }
    }
    
    // Rewire edges with probability beta
    for (auto& edge : edges) {
        if (randProb(gen) < beta) {
            int oldTarget = edge.target;
            std::uniform_int_distribution<int> nodeDist(0, numNodes - 1);
            
            bool validRewire = false;
            int attempts = 0;
            const int maxAttempts = 100;
            
            while (!validRewire && attempts < maxAttempts) {
                int newTarget = nodeDist(gen);
                attempts++;
                
                // Avoid self-loops and duplicates
                if (newTarget != edge.source && 
                    std::find_if(edges.begin(), edges.end(), 
                    [&](const Edge& e) { 
                        return (e.source == edge.source && e.target == newTarget) || 
                               (e.source == newTarget && e.target == edge.source); 
                    }) == edges.end()) {
                    
                    edge.target = newTarget;
                    validRewire = true;
                }
            }
            
            if (!validRewire) {
                edge.target = oldTarget; // Failed to rewire, keep the original
            }
        }
    }
    
    return edges;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <output_file> <graph_type> [parameters]" << std::endl;
        std::cout << "Graph types: random, community, preferential, smallworld" << std::endl;
        std::cout << "Parameters:" << std::endl;
        std::cout << "  random: <numNodes> <numEdges> [minWeight=1.0] [maxWeight=10.0]" << std::endl;
        std::cout << "  community: <numNodes> <numCommunities> <pIntra> <pInter> [minWeight=1.0] [maxWeight=10.0]" << std::endl;
        std::cout << "  preferential: <numNodes> <m0> <m> [minWeight=1.0] [maxWeight=10.0]" << std::endl;
        std::cout << "  smallworld: <numNodes> <k> <beta> [minWeight=1.0] [maxWeight=10.0]" << std::endl;
        return 1;
    }
    
    // Prepend "./graphs/" to output file path
    std::string outputFile = std::string("../inputs/") + argv[1];
    
    std::string graphType = argv[2];
    
    // Set up random number generator with a time-based seed
    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 gen(seed);
    
    std::vector<Edge> edges;
    int numNodes = 0;
    
    try {
        if (graphType == "random") {
            if (argc < 5) {
                throw std::invalid_argument("Random graph requires at least numNodes and numEdges parameters");
            }
            
            numNodes = std::stoi(argv[3]);
            int numEdges = std::stoi(argv[4]);
            double minWeight = (argc > 5) ? std::stod(argv[5]) : 1.0;
            double maxWeight = (argc > 6) ? std::stod(argv[6]) : 10.0;
            
            edges = generateRandomGraph(numNodes, numEdges, minWeight, maxWeight, gen);
        }
        else if (graphType == "community") {
            if (argc < 7) {
                throw std::invalid_argument("Community graph requires numNodes, numCommunities, pIntra, and pInter parameters");
            }
            
            numNodes = std::stoi(argv[3]);
            int numCommunities = std::stoi(argv[4]);
            double pIntra = std::stod(argv[5]);
            double pInter = std::stod(argv[6]);
            double minWeight = (argc > 7) ? std::stod(argv[7]) : 1.0;
            double maxWeight = (argc > 8) ? std::stod(argv[8]) : 10.0;
            
            edges = generateCommunityGraph(numNodes, numCommunities, pIntra, pInter, minWeight, maxWeight, gen);
        }
        else if (graphType == "preferential") {
            if (argc < 6) {
                throw std::invalid_argument("Preferential attachment graph requires numNodes, m0, and m parameters");
            }
            
            numNodes = std::stoi(argv[3]);
            int m0 = std::stoi(argv[4]);
            int m = std::stoi(argv[5]);
            double minWeight = (argc > 6) ? std::stod(argv[6]) : 1.0;
            double maxWeight = (argc > 7) ? std::stod(argv[7]) : 10.0;
            
            edges = generatePreferentialAttachmentGraph(numNodes, m0, m, minWeight, maxWeight, gen);
        }
        else if (graphType == "smallworld") {
            if (argc < 6) {
                throw std::invalid_argument("Small-world graph requires numNodes, k, and beta parameters");
            }
            
            numNodes = std::stoi(argv[3]);
            int k = std::stoi(argv[4]);
            double beta = std::stod(argv[5]);
            double minWeight = (argc > 6) ? std::stod(argv[6]) : 1.0;
            double maxWeight = (argc > 7) ? std::stod(argv[7]) : 10.0;
            
            edges = generateSmallWorldGraph(numNodes, k, beta, minWeight, maxWeight, gen);
        }
        else {
            throw std::invalid_argument("Unknown graph type: " + graphType);
        }
        
        writeGraphToFile(outputFile, numNodes, edges);
        std::cout << "Generated " << graphType << " graph with " << numNodes << " nodes and " 
                  << edges.size() << " edges" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}