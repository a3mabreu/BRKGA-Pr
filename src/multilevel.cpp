#include "csr.hpp"

// Get maximal a Independent Set used in MSH (2001)
std::vector<usize> CSR::maximalIndependentSet() {
    // Vertices are uncolored, colored or forbidden
    std::vector<bool> uncolored(m, true);
    // Colored belongs to the MIS
    std::vector<usize> colored;
    /// Priority queue
    std::vector<VertexCost> gain;
    robin_hood::unordered_map<usize, usize> idx_gain;

    // Prealocating
    colored.reserve(m / 4);
    gain.reserve(m);
    idx_gain.reserve(m);

    // Gain[i] = degree[i]
    for (usize i = 0; i < m; ++i)
        insertHeap(gain, idx_gain, i, degree[i]);

    while (!gain.empty()) {
        const auto [i_max, _] = getFirst(gain, idx_gain);
        uncolored[i_max] = false;
        colored.push_back(i_max);

        // neighbors are removed from the queue (moved into VF)
        for (auto j_idx = row_index[i_max]; j_idx < row_index[i_max + 1]; ++j_idx) {
            const auto j = col_index[j_idx];

            if (uncolored[j]) {  
                removeElement(gain, idx_gain, j);
                uncolored[j] = false;

                // For each j the gain of its uncolored neighbors k are increased by one.
                for (auto k_idx = row_index[j]; k_idx < row_index[j + 1]; ++k_idx) {
                    const auto k = col_index[k_idx];

                    if (uncolored[k]) {
                        usize weight = gain[idx_gain.at(k)].second;
                        ++weight;
                        bubleUp(gain, idx_gain, k, weight);
                    }
                }
            }
        }
    }

    return colored;
}

// Create the Coarse Graph for MSH (2001)
CSR CSR::getCoarseGraph(const std::vector<usize>& mis) {
    // Limits for slow graphs in the new suite of 25 matrices
    constexpr usize MAX_ACC_NEIGHBORS1 = 19'000'000;
    constexpr usize MAX_ACC_NEIGHBORS2 = 47'000'000;
    const usize coarse_m = mis.size();
    const usize max_usize = std::numeric_limits<usize>::max();
    usize nnz = 0;

    // Create mapping from fine vertex index to coarse vertex index
    std::vector<usize> fine_to_coarse(m, max_usize);
    for (usize i = 0; i < coarse_m; ++i)
        fine_to_coarse[mis[i]] = i;

    // The neighbors for each coarse vertex
    std::vector<robin_hood::unordered_set<usize>> coarse_neighbors(coarse_m);

    // Distances works here as a visited time for speed
    distances.assign(m, 0);
    usize current_time = 0;

    // For each MIS vertex, run a BFS (depth < 3) in fine graph
    std::queue<std::pair<usize, usize>> queue;
    for (usize i = 0; i < coarse_m; ++i) {
        const usize start_fine = mis[i];
        // If true, i is not isolated (>1 neighbor)
        bool connected = false;
        ++current_time;

        if (!coarse_neighbors[i].empty())
            connected = true;

        // Store pairs: (current fine vertex, distance)
        queue.push({start_fine, 0});
        distances[start_fine] = current_time;

        while (!queue.empty()) {
            const auto [v, dist] = queue.front();
            queue.pop();

            // If i is connected, the limit can be applyed in slow matrices
            if (connected) {
                const usize n_neighbor = row_index[v + 1] - row_index[v];
                const usize total1 = n_neighbor * (usize)queue.size();
                const usize total2 = n_neighbor * coarse_m;
                if (total1 >= MAX_ACC_NEIGHBORS1 || total2 >= MAX_ACC_NEIGHBORS2) {
                    queue = {};
                    break;
                }
            }

            const usize new_dist = dist + 1;
    
            for (usize idx = row_index[v]; idx < row_index[v + 1]; ++idx) {
                const usize neighbor = col_index[idx];

                // If not visited
                if (distances[neighbor] != current_time) {
                    distances[neighbor] = current_time;
                    // Enqueue only if distance will be <= 3
                    if (new_dist < 3)
                        queue.push({neighbor, new_dist});

                    // If neighbor is in the MIS (and is not the starting vertex)
                    // Add edge from coarse vertex i-th (from start_fine) to the coarse index of 'neighbor'
                    const usize coarse_u = fine_to_coarse[neighbor];
                    if (coarse_u != max_usize && neighbor != start_fine) {
                        coarse_neighbors[i].insert(coarse_u);
                        ++nnz;
                        connected = true;
                        coarse_neighbors[coarse_u].insert(i);
                    }
                }
            }
        }
    }
    
    /// Building the CSR for the coarse graph
    // First, compute row_index for coarse graph
    CSR coarse_csr = CSR(coarse_m, nnz);
    coarse_csr.row_index.resize(coarse_m + 1);
    coarse_csr.labels.resize(coarse_m);
    iota(coarse_csr.labels.begin(), coarse_csr.labels.end(), 0);
    coarse_csr.degree.resize(coarse_m);
    coarse_csr.min_degree = max_usize;
    coarse_csr.max_degree = 0;

    for (usize i = 0; i < coarse_m; ++i) {
        const usize n_size = coarse_neighbors[i].size();
        coarse_csr.degree[i] = n_size;
        coarse_csr.min_degree = std::min(coarse_csr.min_degree, n_size);
        coarse_csr.max_degree = std::max(coarse_csr.max_degree, n_size);

        coarse_csr.row_index[i + 1] = coarse_csr.row_index[i] + n_size;
        std::vector<int> vec(coarse_neighbors[i].begin(), coarse_neighbors[i].end());
        std::sort(vec.begin(), vec.end());
        coarse_csr.col_index.insert(coarse_csr.col_index.end(), vec.begin(), vec.end());
    }

    return coarse_csr;
}

/// Global priority of fine graph
const std::vector<realT> CSR::sloanRefine(const CSR& coarse_csr, const std::vector<usize> mis) {
    // For vertex in the MIS the priority is 
    // the "position" (starting from 1) in coarse after profile reduction
    std::vector<realT> priority(m, 0);
    for (usize i = 0; i < mis.size(); ++i)
        priority[mis[i]] = coarse_csr.labels[i] + 1;

    const robin_hood::unordered_set<usize> mis_set(mis.begin(), mis.end());

    // Others receive the average of the global priority values 
    // of its neighbors that belong to the maximal independent set
    for (usize i = 0; i < m; ++i) {
        if (mis_set.count(i) == 0) {
            realT acc = 0;
            usize msi_neighbors = 0;

            for (usize idx = row_index[i]; idx < row_index[i + 1]; ++idx) {
                const usize j = col_index[idx];

                if (mis_set.count(j) > 0) {
                    acc += priority[j];
                    ++msi_neighbors;
                }
            }

            priority[i] = acc / static_cast<realT>(msi_neighbors);
        }
    }

    return priority;
}

// Multilvel SloanMGPS refinement with SM/MPG in the coarsest graph
// ALGO_BASE:  0 = SM; 1 = MPG
void CSR::msW(const std::vector<realT>& priority, const usize ALGO_BASE) {
    constexpr float MAX_RATIO = 0.8;
    thread_local usize level = 0;
    constexpr usize MAX_LEVEL = 1;

    // Base case: if the graph is sufficiently small
    if (level >= MAX_LEVEL || m <= 2) {
        // SloanMGPS or MPG in the coarsest graph
        if (priority.empty()) {
            if (ALGO_BASE == 0)
                enhancedSloanMGPSPriority({});
            else 
                mpg();
        } else {
            enhancedSloanMGPSPriority(priority);
        }

        return;
    }

    // Coarsening phase: obtain the coarse graph via MIS selection
    const std::vector<usize> mis = maximalIndependentSet();
    CSR coarse_csr = getCoarseGraph(mis);

    // Check ratio
    const float ratio = static_cast<float>(coarse_csr.m) / static_cast<float>(m);
    if (ratio > MAX_RATIO) {
        // SloanMGPS or MPG in the coarsest graph
        if (priority.empty()) {
            if (ALGO_BASE == 0)
                enhancedSloanMGPSPriority({});
            else 
                mpg();
        } else {
            enhancedSloanMGPSPriority(priority);
        }

        return;
    }

    /// First descent: recursive call down to the coarser level.
    ++level;
    if (priority.empty()) {
        coarse_csr.msW({}, ALGO_BASE);
    } else {
        coarse_csr.msW(priority, ALGO_BASE);
    }
    --level;

    // First prolongation and refinement step
    const std::vector<realT>& priority1 = sloanRefine(coarse_csr, mis);
    enhancedSloanMGPSPriority(priority1);

    /// The W additional re-coarsening of the refined graph
    // No need to re-compute MIS (it's the same)
    static std::vector<realT> pri(m); // Initialize with the correct size
    std::transform(labels.begin(), labels.end(), pri.begin(), [](usize val) { return static_cast<realT>(val); });
    
    ++level;
    coarse_csr.msW(pri, ALGO_BASE);
    --level;

    // Second prolongation/refinement step
    const std::vector<realT>& priority2 = sloanRefine(coarse_csr, mis);
    enhancedSloanMGPSPriority(priority2);
}

void CSR::msWConstrutiveSM(const std::vector<realT>& priority) {
    static constexpr float MAX_RATIO = 0.8;
    thread_local usize level = 0;
    static const usize MAX_LEVEL = 1;

    // Base case: if the graph is sufficiently small
    if (level >= MAX_LEVEL || m <= 2) {
        if (priority.empty()) {
            constructiveNSloanMGPS(realZeroOneInclusive());
        } else {
            enhancedSloanMGPSPriority2(priority);
        }

        return;
    }

    // Coarsening phase: obtain the coarse graph via MIS selection
    const std::vector<usize> mis = maximalIndependentSet();
    CSR coarse_csr = getCoarseGraph(mis);

    // Check ratio
    const float ratio = static_cast<float>(coarse_csr.m) / static_cast<float>(m);
    if (ratio > MAX_RATIO) {
        if (priority.empty()) {
            constructiveNSloanMGPS(realZeroOneInclusive());
        } else {
            enhancedSloanMGPSPriority2(priority);
        }

        return;
    }

    /// First descent: recursive call down to the coarser level.
    ++level;
    if (priority.empty()) {
        coarse_csr.msWConstrutiveSM({});
    } else {
        coarse_csr.msWConstrutiveSM(priority);
    }
    --level;

    // First prolongation and refinement step
    const std::vector<realT>& priority1 = sloanRefine(coarse_csr, mis);
    enhancedSloanMGPSPriority2(priority1);

    /// The W additional re-coarsening of the refined graph
    // No need to re-compute MIS (it's the same)
    static std::vector<realT> pri(m); // Initialize with the correct size
    std::transform(labels.begin(), labels.end(), pri.begin(), [](usize val) { return static_cast<realT>(val); });
    
    ++level;
    coarse_csr.msWConstrutiveSM(pri);
    --level;

    // Second prolongation/refinement step
    const std::vector<realT>& priority2 = sloanRefine(coarse_csr, mis);
    enhancedSloanMGPSPriority2(priority2);
}
