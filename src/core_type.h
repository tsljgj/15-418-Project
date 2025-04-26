// core_type.h
#ifndef CORE_TYPE_H
#define CORE_TYPE_H

// Enum for specifying CPU core type
enum CoreType {
    ANY_CORE,  // No specific core preference
    P_CORE,    // Performance core
    E_CORE     // Efficiency core
};

// Set thread affinity to a specific CPU ID
void setThreadAffinityToCpu(int cpuId);

// Get the appropriate CPU ID range for a core type
std::vector<int> getCpuIdsForCoreType(CoreType coreType);

// Set thread affinity based on core type
void setThreadAffinityCoreType(CoreType coreType);

// Helper function to set thread affinity to specific core type (legacy function)
void setCoreAffinity(CoreType coreType);

#endif // CORE_TYPE_H