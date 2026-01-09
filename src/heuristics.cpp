#include "csr.hpp"

// Algorithm Sloan-MGPS(1999) for profile reduction
// For original Sloan (1989) use (w1=2,w2=1)
void CSR::sloanMGPS(realT w1, realT w2, const bool normalized) {
    // Data of vertices
    std::vector<VertexDataReal> vertices(m);
    /// Priority queue (vertices active or preactive)
    std::vector<VertexCostReal> Q;
    robin_hood::unordered_map<usize, usize> idx_Q;

    // Misc
    usize current_label = 0;
    realT max_real = std::numeric_limits<realT>::max();

    /// Step 1 and 2 - Pseudo-peripheral nodes s and e
    auto [s, e] = mgpsPP();

    // Compute distances from e (in vector *distance*)
    bfs(e);

    /// Weights
    // Here I'm using the weights as in Reid and Scott (99) and Hu and Scott (2001)
    if (normalized) {
        // max_d is distance from s to e (undirected graph!)
        const realT max_d = static_cast<realT>(distances[s]);
        realT norm = max_d / static_cast<realT>(max_degree);
        // If the pseudo-diameter < max_degree, norm is one
        norm = std::max(static_cast<realT>(1.0), norm);
        w1 *= norm;
    }

    /// Step 3
    for (usize i = 0; i < m; ++i) {
        vertices[i].priority = -w1 * (degree[i] + 1) + w2 * distances[i] ;
        vertices[i].cur_degree = degree[i];
        vertices[i].status = VertexStatus::Inactive;
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
            // The increment in priority means that current degree is reduced by 1
            for (usize j_idx = row_index[i]; j_idx < row_index[i + 1]; ++j_idx) {
                const usize j = col_index[j_idx];

                // vertices[j].priority += w1;
                --vertices[j].cur_degree;
                if (vertices[j].cur_degree > 0) {
                    vertices[j].priority += w1;
                } else {
                    vertices[j].priority = max_real;
                }

                if (vertices[j].status == VertexStatus::Inactive) {
                    vertices[j].status = VertexStatus::Preactive;
                    insertHeapReal(Q, idx_Q, j, vertices[j].priority);
                } else if (vertices[j].status != VertexStatus::Postactive) {
                    // Update priority in Q
                    // if (idx_Q.count(j) > 0) // Shouldn't be necessary
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
                
                // vertices[j].priority += w1;
                --vertices[j].cur_degree;
                if (vertices[j].cur_degree > 0) {
                    vertices[j].priority += w1;
                } else {
                    vertices[j].priority = max_real;
                }

                bubleUpReal(Q, idx_Q, j, vertices[j].priority);

                // For each k in Adj(j)
                for (usize k_idx = row_index[j]; k_idx < row_index[j + 1]; ++k_idx) {
                    const usize k = col_index[k_idx];

                    if (vertices[k].status != VertexStatus::Postactive) {
                        // vertices[k].priority += w1;
                        --vertices[k].cur_degree;
                        if (vertices[k].cur_degree > 0) {
                            vertices[k].priority += w1;
                        } else {
                            vertices[k].priority = max_real;
                        }

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

// SloanMGPS (MC60) with global priority function (4) of the hibryd Sloan 
// from Reid and Scott (1999) and (2.5) of Hu and Scoot (2001)
void CSR::sloanMGPSPriority(const realT w1, const realT w2, const std::vector<realT>& priority) {
    // Data of vertices
    std::vector<VertexDataReal> vertices(m);
    /// Priority queue (vertices active or preactive)
    std::vector<VertexCostReal> Q;
    robin_hood::unordered_map<usize, usize> idx_Q;
    // Misc
    usize current_label = 0;
    realT max_real = std::numeric_limits<realT>::max();

    /// Step 1 and 2 - Pseudo-peripheral nodes s and e
    auto [s, e] = mgpsPP();

    // Compute distances from e (in vector *distance*)
    bfs(e);

    // h is the level-set depth (distance from e to s)
    const realT h = distances[s];

    /// The priority is for the MSH method
    // Here I'm using the weights as in Reid and Scott (99) and Hu and Scott (2001)
    // P_i = -W1 * c_i -W2 (h/n) p_i
    const realT nu = w2 * (h / static_cast<realT>(m));
    if (!priority.empty()) {
        for (usize i = 0; i < m; ++i) {
            vertices[i].priority = -w1 * (degree[i] + 1) -nu * priority[i];
            vertices[i].cur_degree = degree[i];
            vertices[i].status = VertexStatus::Inactive;
        }
    } else {
        for (usize i = 0; i < m; ++i) {
            vertices[i].priority = -w1 * (degree[i] + 1) -nu * distances[i];
            vertices[i].cur_degree = degree[i];
            vertices[i].status = VertexStatus::Inactive;
        }
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
            // The increment in priority means that current degree is reduced by 1
            for (usize j_idx = row_index[i]; j_idx < row_index[i + 1]; ++j_idx) {
                const usize j = col_index[j_idx];

                // vertices[j].priority += w1;
                --vertices[j].cur_degree;
                if (vertices[j].cur_degree > 0) {
                    vertices[j].priority += w1;
                } else {
                    vertices[j].priority = max_real;
                }

                if (vertices[j].status == VertexStatus::Inactive) {
                    vertices[j].status = VertexStatus::Preactive;
                    insertHeapReal(Q, idx_Q, j, vertices[j].priority);
                } else if (vertices[j].status != VertexStatus::Postactive) {
                    // Update priority in Q
                    // if (idx_Q.count(j) > 0) // Shouldn't be necessary
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
                
                // vertices[j].priority += w1;
                --vertices[j].cur_degree;
                if (vertices[j].cur_degree > 0) {
                    vertices[j].priority += w1;
                } else {
                    vertices[j].priority = max_real;
                }

                bubleUpReal(Q, idx_Q, j, vertices[j].priority);

                // For each k in Adj(j)
                for (usize k_idx = row_index[j]; k_idx < row_index[j + 1]; ++k_idx) {
                    const usize k = col_index[k_idx];

                    if (vertices[k].status != VertexStatus::Postactive) {
                        // vertices[k].priority += w1;
                        --vertices[k].cur_degree;
                        if (vertices[k].cur_degree > 0) {
                            vertices[k].priority += w1;
                        } else {
                            vertices[k].priority = max_real;
                        }

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


// Enhanced Sloan(1993) with MGPS (MC60)
// This method uses the priority function (4) of the enhanced Sloan from Reid and Scott (1999) and (2.5) of Hu and Scoot (2001) with 2 pairs of weights [(2,1) and (16,1)]
void CSR::enhancedSloanMGPS() {
    // Considers natural labeling
    std::vector<usize> best_labels = labels;
    evaluateProfile();
    unsigned long best_profile = profile;

    // Try (2,1) and (16,1) for the SloanMGPS
    sloanMGPS(static_cast<realT>(2), static_cast<realT>(1), false);
    evaluateProfile();
    // Update the best
    if (profile < best_profile) {
        best_profile = profile ;
        std::swap(labels, best_labels);
    }

    sloanMGPS(static_cast<realT>(16), static_cast<realT>(1), false);
    evaluateProfile();
    // If do not improves, revert to the best
    if (profile > best_profile) {
        profile = best_profile;
        std::swap(labels, best_labels);
    }
}

// Enhanced Sloan(1993) with MGPS (MC60)
// This method uses the priority function (4) of the enhanced Sloan from Reid and Scott (1999) and (2.5) of Hu and Scoot (2001) with 2 pairs of weights [(2,1) and (16,1)]
void CSR::enhancedSloanMGPSPriority(const std::vector<realT>& priority) {
    // Considers natural labeling
    std::vector<usize> best_labels = labels;
    evaluateProfile();
    unsigned long best_profile = profile;

    sloanMGPSPriority(static_cast<realT>(2), static_cast<realT>(1), priority);
    evaluateProfile();
    // Update the best
    if (profile < best_profile) {
        best_profile = profile ;
        std::swap(labels, best_labels);
    }

    sloanMGPSPriority(static_cast<realT>(16), static_cast<realT>(1), priority);
    evaluateProfile();
    // If do not improves, revert to the best
    if (profile > best_profile) {
        profile = best_profile;
        std::swap(labels, best_labels);
    } 
}

// For populational algorithms (no natual labeling for diversity)
void CSR::enhancedSloanMGPSPriority2(const std::vector<realT>& priority) {
    sloanMGPSPriority(static_cast<realT>(2), static_cast<realT>(1), priority);
    evaluateProfile();
    unsigned long best_profile = profile;
    std::vector<usize> best_labels = labels;

    sloanMGPSPriority(static_cast<realT>(16), static_cast<realT>(1), priority);
    evaluateProfile();
    // If do not improves, revert to the best
    if (profile > best_profile) {
        profile = best_profile;
        std::swap(labels, best_labels);
    } 
}

// Algorithm MPG(1993) for profile reduction
void CSR::mpg() {
    // (a) l is the list of nodes which have already been labelled, i.e. lj Є l is the node with label j.
    // I will use labels and check if is not infinity
    const usize max_usize = std::numeric_limits<usize>::max();
    labels.assign(m, max_usize);
    // (d) The set Y contains the elements of l or q
    // (e) The current degree d of a node at a given stage of labelling is equal to the number of nodes still not labelled adjacent to i but not belonging to q
    // d[i] = |Adj(i) - (Adj(i) \cap Y)|
    std::vector<usize> d;
    // (f) The priority of a node i at a given stage of labelling for a end node e is defined by: 
    // p[i] = d_e[i] - 2 * d[i] 
    std::vector<int> p(m);
    // Number of connections to q
    std::vector<usize> a;
    // Misc
    std::vector<VertexCost> aux_T;
    robin_hood::unordered_map<usize, usize> aux_idx_T;
    usize current_label = 0;
    usize N; // used in step 4 and 6


    ///////////
    /// (1) (Entry) Enter with the pseudo-peripheral nodes s, e and the rooted level structure
    auto [s, e] = sloanPP();
    bfs(e);
    // Distances d_e of nodes from the vertex e
    std::vector<usize> d_e = std::move(distances);

    /// (2) (Initialize priority queues) Set q = Ø and t = Ø
    robin_hood::unordered_map<usize, usize> idx_T;
    robin_hood::unordered_map<usize, usize> idx_Q;
    // Container of pairs (<vertex, priority>) for the heap.
    // (c) q is the priority queue (max heap) of nodes that belong to t or are eligible to be inserted in t.
    std::vector<VertexCost> Q;
    // T is the priority queue (max heap) of nodes which are eligible to be labelled.
    std::vector<VertexCost> T;


    /// (3) (Initialize current degrees d, priorities and numbers of connections to q
    d = degree;
    for (usize i = 0; i < m; ++i)
        p[i] = d_e[i] - 2 * d[i];
    a.assign(m, 0);
    // The first vertex to be included in Q is s
    N = s;

    do {
        // (4) (Insert a new node in q)  If |q| = 0, then set q = {s}
        // Else, insert at the end of Q a node n \in Adj(t) - (Adj(t) \cap Y) such that this formula is maximum:
        // π_n = 2p_n + 2p_max (Adj({n}) \cap t) + 3a+n
        int pi_max = std::numeric_limits<int>::min();

        for (const VertexCost& v : T) {
            const usize u = v.first;
            for (auto j_idx = row_index[u]; j_idx < row_index[u + 1]; ++j_idx) {
                const auto n = col_index[j_idx];

                // Remove elements in Y (elements in L or in Q)
                if (labels[n] != max_usize || idx_Q.count(n) > 0)
                    continue;

                // Get p_max of Adj of n in T
                int p_max = std::numeric_limits<int>::min();
                for (auto j_idx = row_index[n]; j_idx < row_index[n + 1]; ++j_idx) {
                    const auto adj_n = col_index[j_idx];

                    if (idx_T.count(adj_n) > 0 && p[adj_n] > p_max)
                        p_max = p[adj_n];
                }

                const int pi = 2 * p[n] + 2 * p_max + 3 * a[n];
                // Get maximum priority
                if (pi > pi_max) {
                    pi_max = pi;
                    N = n;
                }
            }                
        }

        // Insert max priority element in Q not in Y
        if (labels[N] == max_usize && idx_Q.count(N) == 0) {
            insertHeap(Q, idx_Q, N, p[N]);

            /// (6) (Insert nodes in T) Scan Adj({N}) ∩ q and insert at the end of T nodes for which d[i] = 1 and NOT IN T.
            for (auto j_idx = row_index[N]; j_idx < row_index[N + 1]; ++j_idx) {
                const auto j = col_index[j_idx];
                // Updade d and a for Adjs of N (inserted in Q) (Step 5)
                --d[j];
                ++a[j];
                // Update priority (Step 5)
                p[j] = d_e[j] - 2 * d[j];
                // Restore heap property using bubleUp
                if (idx_Q.count(j) > 0)
                    bubleUp(Q, idx_Q, j, p[j]);
                if (idx_T.count(j) > 0) {
                    bubleUp(T, idx_T, j, p[j]);
                } else {
                    // if j is in Q, not in T and d[j] == 1
                    if (idx_Q.count(j) > 0 && d[j] == 1) {
                        // Insert in heap, update idx and heapify
                        insertHeap(T, idx_T, j, p[j]);
                    }
                }
                
            }
        }

        /// (7) (Label nodes of T) Scan t in ascending order extracting from T, deleting in Q, and inserting
        // at the end of L the nodes iЄt for which d[i] = 0.
        aux_T.clear();
        aux_idx_T.clear();
        while (!T.empty()) {
            /// Get first, do min-heapify and update idx_H 
            const auto& v = getFirst(T, idx_T);
            const usize i = v.first;
            const usize p_i = v.second;

            // Label the vertex if d[i] == 0
            if (d[i] == 0) {
                labels[i] = current_label;
                ++current_label;

                // Remove from Q
                if (idx_Q.count(i) > 0) {
                    removeElement(Q, idx_Q, i);
                    // Decrease connections to Q for Adj(i)
                    for (auto j_idx = row_index[i]; j_idx < row_index[i + 1]; ++j_idx) {
                        const auto j = col_index[j_idx];
                        --a[j];
                    }
                }
            } else {
                // Must be inserted in T again if d[i] <=1
                if (d[i] <= 1)
                    insertHeap(aux_T, aux_idx_T, i, p_i);
            }
        }

        /// Restore T
        if (!aux_T.empty()) {
            T = std::move(aux_T);
            idx_T = std::move(aux_idx_T);
        }

        /// (8) (Build new T)
        // If |T| = 0, scan Q in ascending order and insert at the end of T nodes iЄq for which p[i] >= p_max(Q) — 1.
        if (T.empty() && !Q.empty()) {
            // The top is p_max(Q)
            int p_max_q = (Q[0].second) - 1;
            // Insert in T
            for (const auto &v : Q) {
                const usize i = v.first;
                const int pri = v.second;

                if (pri >= p_max_q) {
                    insertHeap(T, idx_T, i, pri);
                }
            }
        }
    } while (current_label < m);
}
