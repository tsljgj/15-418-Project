// core_type.cpp
#include "core_type.h"
#include <iostream>
#include <thread>
#include <vector>
#include <windows.h>
#include <omp.h>
#include <chrono>

// Set thread affinity to a specific CPU ID
void setThreadAffinityToCpu(int cpuId) {
    HANDLE currentThread = GetCurrentThread();
    DWORD_PTR affinityMask = (static_cast<DWORD_PTR>(1) << cpuId);
    DWORD threadId = GetCurrentThreadId();
    
    std::cout << "Setting thread " << threadId << " to CPU " << cpuId << std::endl;
    
    BOOL result = SetThreadAffinityMask(currentThread, affinityMask);
    
    if (result == 0) {
        std::cerr << "Error setting thread affinity to CPU " << cpuId 
                  << ": " << GetLastError() << " for Thread ID: " << threadId << std::endl;
    } else {
        std::cout << "Thread " << threadId << " successfully assigned to CPU " << cpuId << std::endl;
                //   << " (SetThreadAffinityMask took " << elapsed_seconds.count() << " seconds)" << std::endl;
    }
}

// Get the appropriate CPU ID range for a core type
std::vector<int> getCpuIdsForCoreType(CoreType coreType) {
    std::vector<int> cpuIds;
    int numCPUs = std::thread::hardware_concurrency();
    
    switch (coreType) {
        case P_CORE:
            // P-cores are 0-15 on your system (8 P cores, two threads per P core)
            for (int i = 0; i < 16 && i < numCPUs; i++) {
                cpuIds.push_back(i);
            }
            break;
            
        case E_CORE:
            // E-cores are 16-31 on your system
            for (int i = 16; i < 32 && i < numCPUs; i++) {
                cpuIds.push_back(i);
            }
            break;
            
        case ANY_CORE:
        default:
            // All available cores
            for (int i = 0; i < numCPUs; i++) {
                cpuIds.push_back(i);
            }
            break;
    }
    
    return cpuIds;
}

// Set thread affinity based on core type
void setThreadAffinityCoreType(CoreType coreType) {
    // If ANY_CORE, don't set any specific affinity
    if (coreType == ANY_CORE) {
        std::cout << "No specific core affinity set - using any available core (system decides)\n";
        return;
    }
    
    std::vector<int> cpuIds = getCpuIdsForCoreType(coreType);
    if (cpuIds.empty()) {
        std::cerr << "No CPUs found for the specified core type\n";
        return;
    }
    
    HANDLE currentThread = GetCurrentThread();
    DWORD_PTR affinityMask = 0;
    
    // Create affinity mask for all CPUs of the requested type
    for (int cpuId : cpuIds) {
        affinityMask |= (static_cast<DWORD_PTR>(1) << cpuId);
    }
    
    DWORD threadId = GetCurrentThreadId();
    std::cout << "Setting thread " << threadId << " to ";
    if (coreType == P_CORE) {
        std::cout << "P-cores (CPUs 0-15)\n";
    } else if (coreType == E_CORE) {
        std::cout << "E-cores (CPUs 16-31)\n";
    }
    
    BOOL result = SetThreadAffinityMask(currentThread, affinityMask);
    if (result == 0) {
        std::cerr << "Error setting thread affinity: " << GetLastError() 
                 << " for Thread ID: " << threadId << std::endl;
    } else {
        std::cout << "Affinity set successfully for Thread ID: " << threadId << std::endl;
        
        // Get and print the actual affinity mask that was set
        DWORD_PTR processAffinityMask, systemAffinityMask;
        if (GetProcessAffinityMask(GetCurrentProcess(), &processAffinityMask, &systemAffinityMask)) {
            std::cout << "Process affinity mask: 0x" << std::hex << processAffinityMask << std::dec << std::endl;
        }
    }
}

// Legacy function to maintain compatibility with existing code
void setCoreAffinity(CoreType coreType) {
    // Simply call the new function
    setThreadAffinityCoreType(coreType);
}

// Assign specific P-cores and E-cores for parallel execution
std::vector<int> assignParallelCores(int pCoreCount, int eCoreCount) {
    std::vector<int> coreAssignments;
    
    // Validate core counts
    if (pCoreCount > 8) {
        std::cerr << "Error: Requested P-core count (" << pCoreCount << ") exceeds maximum (8)\n";
        return coreAssignments;
    }
    
    if (eCoreCount > 16) {
        std::cerr << "Error: Requested E-core count (" << eCoreCount << ") exceeds maximum (16)\n";
        return coreAssignments;
    }
    
    // Assign P-cores (use odd indices to avoid hyperthreading)
    for (int i = 0; i < pCoreCount; i++) {
        coreAssignments.push_back(i * 2 + 1); // 1, 3, ...
    }
    
    // Assign E-cores (16-31)
    for (int i = 0; i < eCoreCount; i++) {
        coreAssignments.push_back(16 + i); // 16, 17, 18, ...
    }
    
    return coreAssignments;
}