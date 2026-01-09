#ifndef MY_TYPES_H
#define MY_TYPES_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <queue>
#include <stack>
#include <climits>
#include <numeric> // for std::iota
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <chrono>
#include <variant>
#include "robin_hood.h"

#define usize unsigned int
#define realT double
using ParamValue = std::variant<int, std::string, realT, int64_t>;

// Status of vertices used in Sloan algorithm
enum class VertexStatus {
    Postactive, // The vertex is already numbered
    Active,     // Unumbered vertex that are adjacent to some numbered vertex
    Preactive,  // Adjacent to active vertex but are neither active nor numbered
    Inactive    // None of the above
};

// Information of vertices used in Sloan algorithm
struct VertexData {
  VertexStatus status;
  int priority;
};

struct VertexDataReal {
  realT priority;
  usize cur_degree;
  VertexStatus status;
};

// pair<vertex, priority> (used in the max heap)
typedef std::pair<usize, int> VertexCost;
typedef std::pair<usize, realT> VertexCostReal;
struct CompareCost {
    bool operator()(const VertexCost& a, const VertexCost& b) {
        return a.second < b.second;
    }
};
struct CompareCostReal {
    bool operator()(const VertexCost& a, const VertexCost& b) {
        return a.second < b.second;
    }
};

using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

using Duration = int64_t;

// Solution Struct for BRKGA
struct SolutionRK {
    unsigned long profile; // Solution objective function value
    std::vector<usize> labels; // labels == Current solution
    std::vector<realT> random_keys; // Representation using Random Keys

    bool operator<(const SolutionRK& other) const {
    if (profile < other.profile)
        return true;
    else 
        return false;
    }
};

struct SolutionDR {
    unsigned long profile;
    std::vector<usize> labels;
    std::vector<usize> rho;
    // For sorting and comparison
    bool operator<(const SolutionDR& other) const {
        return profile < other.profile;
    }
    bool operator>(const SolutionDR& other) const {
        return profile > other.profile;
    }
    bool operator==(const SolutionDR& other) const {
        return profile == other.profile;
    }
};

// Structure to hold element value and original index
struct IndexedElement {
    realT value;
    usize index;
};

struct Element {
    usize i;
    usize j;
    realT v;
};

#endif