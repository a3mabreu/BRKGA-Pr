#ifndef CSR_H
#define CSR_H

#include "types.hpp"
#include "misc.hpp"
#include "max_heap_robin_hood.hpp"
#include "max_heap_real.hpp"

class CSR {
public:
    std::vector<IndexedElement> indexed_rk; // Store the RK and its original index
    std::vector<realT> random_keys; // Representation using Random Keys
    std::vector<realT> tmp_rk;  // Aux vector used in the encoder (random_keys.cpp) and in psiVertices
    /*************/
    std::vector<usize> col_index; // Column indices of non-zero entries
    std::vector<usize> row_index; // Row index pointers
    std::vector<usize> labels; // Labels of vertices
    std::vector<usize> degree; // Degre of vertices
    std::vector<char> visited; // Visited vertices
    std::vector<usize> distances; // Distances
    std::vector<usize> reduced_n; // Reduced neighbourhood for LS
    // Vertices affected in swaping and updating profile 
    TimePoint t_start; // Start time
    std::chrono::seconds max_time; // Maximum execution time
    unsigned long profile; // Profile
    unsigned long best_profile = std::numeric_limits<unsigned long>::max(); // Best Profile so far
    realT alpha; // alpha for the msW with construtctiveMPG
    usize m; // Number of rows
    usize n_nz; // Number of non-zero entries
    usize max_degree = 0; // Maximum degree of instance
    usize min_degree; // Maximum degree of instance
    bool symmetric;

    // Simple construtor
    CSR(usize rows, usize nnz);
    // Constructor overloaded for reading .mtx files
    CSR(const std::string& path, bool f_symmetric = false);

    // Evaluate profile
    void evaluateProfile();

    /// Pseudoperipheral vertex
    // Get a pair of pseudoperipheral (s, e) (Sloan's algorithm)
    std::pair<usize, usize> sloanPP();
    // Get a pair of pseudoperipheral (start, end) MGPS algorithm)
    std::pair<usize, usize> mgpsPP();

    /// Profile reduction heuristics
    // Algorithm MPG(1993) for profile reduction
    void mpg();
    // Algorithm Sloan-MGPS(1999) for profile reduction
    void sloanMGPS(realT w1 = 2, realT w2 = 1, const bool normalized = false);
    // Enhanced version with 2 pairs of weights
    void enhancedSloanMGPSPriority(const std::vector<realT>& priority = {});
    // For populational algorithms (no natual labeling for diversity)
    void enhancedSloanMGPSPriority2(const std::vector<realT>& priority = {});
    // This method uses the priority function (4) of the enhanced Sloan from Reid and Scott (1999) and (2.5) of Hu and Scoot (2001)
    void enhancedSloanMGPS();

    /// Constructive methods
    // Randomized semi-greed constructive method based on SloanMGPS
    void constructiveNSloanMGPS(const realT ALPHA);
    void msWConstrutiveSM(const std::vector<realT>& priority = {});

    /// Multilevel approaches
    // Multilvel SloanMGPS refinement with SM/MPG in the coarsest graph
    void msW(const std::vector<realT>& priority = {}, const usize ALGO_BASE = 0);
    /// Labeling the fine graph (MSH)
    const std::vector<realT> sloanRefine(const CSR& coarse_csr, const std::vector<usize> mis);
    // Get maximal a Independent Set used in MSH (2001)
    std::vector<usize> maximalIndependentSet();
    // Create the Coarse Graph for MSH (2001)
    CSR getCoarseGraph(const std::vector<usize>& mis);
    // SloanMGPS (MC60) with global priority function
    void sloanMGPSPriority(const realT w1 = 2, const realT w2 = 1, const std::vector<realT>& priority = {});

    /// Misc
    // Get vertices from the last level structure and eccentricity
    std::pair<std::vector<usize>, usize> getLastLevelAndEccentricity(usize v);
    // Get eccentricity and level width (max) of v
    std::pair<usize, usize> getEccentricityNWidth(usize v);
    // Get the diameter of the graph
    usize getDiameter();
    // Recursive Breadth-First Search (DFS) 
    void bfs(usize v);
    // Verify if solution is viable (LABELS only)
    bool isFeasible() const;
};

// https://math.nist.gov/MatrixMarket/formats.html
#endif