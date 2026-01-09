#include "csr.hpp"

// Randomized constructive method based on SloanMGPS
// !!! ALPHA must be in the [0, 1] interval !!!
void CSR::constructiveNSloanMGPS(const realT ALPHA) {
    // Data of vertices
    std::vector<VertexDataReal> vertices(m);
    /// Priority queue (vertices active or preactive)
    std::vector<VertexCostReal> Q;
    robin_hood::unordered_map<usize, usize> idx_Q;
    // Misc
    usize current_label = 0;

    /// Step 1 and 2 - Pseudo-peripheral nodes s and e
    auto [s, e] = mgpsPP();
    // Compute distances from e (in vector *distance*)
    bfs(e);

    // Weights
    // max_d is distance from s to e (undirected graph!)
    const realT max_d = static_cast<realT>(distances[s]);
    realT norm = max_d / static_cast<realT>(max_degree);
    // If the pseudo-diameter is less than the max_degree, norm is one
    norm = std::max(static_cast<realT>(1.0), norm);
    const realT w1 = ALPHA;
    const realT w2 = norm * (static_cast<realT>(1) - w1);

    /// Step 3
    for (usize i = 0; i < m; ++i) {
        vertices[i].status = VertexStatus::Inactive;
        vertices[i].priority = w1 * distances[i] - w2 * (degree[i] + 1);
    }
    
    /// Step 4
    insertHeapReal(Q, idx_Q, s, vertices[s].priority);
    vertices[s].status = VertexStatus::Preactive;

    /// Step 5
    while (!Q.empty()) {
        /// Steps 6 and 7
        const auto [i, p_i] = getFirstReal(Q, idx_Q);

        // If node i is preactive
        if (vertices[i].status == VertexStatus::Preactive) {
            // For each j in Adj(i) increment priority
            for (usize j_idx = row_index[i]; j_idx < row_index[i + 1]; ++j_idx) {
                const usize j = col_index[j_idx];

                vertices[j].priority += w2;

                if (vertices[j].status == VertexStatus::Inactive) {
                    vertices[j].status = VertexStatus::Preactive;
                    insertHeapReal(Q, idx_Q, j, vertices[j].priority);
                } else if (vertices[j].status != VertexStatus::Postactive) {
                    // Update priority in Q
                    bubleUpReal(Q, idx_Q, j, vertices[j].priority);
                }
            }
        }

        /// Step 8 - Label and set as postactive
        labels[i] = current_label;
        ++current_label;
        vertices[i].status = VertexStatus::Postactive;

        /// Step 9 - (Update priorities and queue) Examine each node j which is adjacent to node i.
        for (usize j_idx = row_index[i]; j_idx < row_index[i + 1]; ++j_idx) {
            const usize j = col_index[j_idx];

            if (vertices[j].status == VertexStatus::Preactive) {
                vertices[j].status = VertexStatus::Active;
                vertices[j].priority += w2;

                bubleUpReal(Q, idx_Q, j, vertices[j].priority);

                // For each k in Adj(j)
                for (usize k_idx = row_index[j]; k_idx < row_index[j + 1]; ++k_idx) {
                    const usize k = col_index[k_idx];

                    if (vertices[k].status != VertexStatus::Postactive) {
                        vertices[k].priority += w2;

                        if (vertices[k].status == VertexStatus::Inactive) {
                            vertices[k].status = VertexStatus::Preactive;
                            insertHeapReal(Q, idx_Q, k, vertices[k].priority);
                        } else { 
                            // Should be preactive or active
                            bubleUpReal(Q, idx_Q, k, vertices[k].priority);
                        }
                    }
                }
            }
        }
    }
}
