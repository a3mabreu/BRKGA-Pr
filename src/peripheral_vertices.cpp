#include "csr.hpp"

// Get a pair of pseudoperipheral (start, end) for MPG (Sloan's algorithm)
std::pair<usize, usize> CSR::sloanPP() {
    // List of nodes in last level
    std::vector<std::pair<usize, usize>> f;
    // Width of the level structure rooted at e
    usize width_e;
    // pair os vertices starting and end vertices
    usize s, e = 0;
    bool end = false; // termination flag

    /// Select a random s with the smallest degree
    std::vector<usize> candidates;
    for (auto it = degree.begin(); it != degree.end(); ++it) {
        if (*it == min_degree) {
            s = std::distance(degree.begin(), it);
            candidates.push_back(s);
        }
    }
    const usize x = usizeRandomNumber(0, candidates.size() - 1);
    s = candidates[x];

    while (!end) {
        // Generate rooted level structure rooted at s
        auto [last_level_s, eccentricity_s] = getLastLevelAndEccentricity(s);
        end = true;
        /// (Step 3) Sort the last level in ascending sequence of their degrees
        f.clear();
        for (usize i = 0; i < last_level_s.size(); ++i) {
            const usize v = last_level_s[i];
            f.push_back({v, degree[v]});
        }
        // Sort by degree (second element)
        std::sort(f.begin(), f.end(), [](const std::pair<usize, usize> &a, const std::pair<usize, usize> &b) {
            return a.second < b.second;
        });

        /// (Shrink the last level) Scan the sorted level and form a list of nodes f containing only one node of each degree
        // Compare degrees for uniqueness (it_end_f is where f ends)
        const auto it_end_f = std::unique(f.begin(), f.end(), [](const std::pair<usize, usize> &a, const std::pair<usize, usize> &b) {
            return a.second == b.second;
        });

        /// Set w(e) = infinity
        width_e = std::numeric_limits<usize>::max();

        /// (Test for termination) For each node i in f in order of ascending degree generate root level
        for (auto it = f.begin(); it != it_end_f; ++it) {
            const usize i = (*it).first;
            auto [eccentricity_i, width_i] = getEccentricityNWidth(i);

            // If h(i) > h(s) and w(i) < w(e) then set s = i and go to step (3)
            // Else, if w(i) < w(e) then set e = i and w(e) = w(i)
            if (eccentricity_i > eccentricity_s && width_i < width_e) {
                s = i;
                end = false;
                break; // Go to step 3
            } else if (width_i < width_e) {
                e = i;
                width_e = width_i;
            }

        }
    }

    return std::make_pair(s, e);
}


// Get a pair of pseudoperipheral (start, end) MGPS algorithm)
std::pair<usize, usize> CSR::mgpsPP() {
    // List of nodes in last level
    std::vector<usize> last_level;

    // Considered vertices in last level
    std::unordered_set<usize> considered_vertices;
    
    // Width of the level structure rooted at e
    usize width_e;
    // pair os vertices starting and end vertices
    usize s, e = 0;
    bool end = false; // termination flag

    /// Select a random s with the smallest degree
    std::vector<usize> candidates;
    for (auto it = degree.begin(); it != degree.end(); ++it) {
        if (*it == min_degree) {
            s = std::distance(degree.begin(), it);
            candidates.push_back(s);
        }
    }
    const usize x = usizeRandomNumber(0, candidates.size() - 1);
    s = candidates[x];

    while (!end) {
        end = true;
        // Generate rooted level structure rooted at s
        auto [last_level_s, eccentricity_s] = getLastLevelAndEccentricity(s);

        // Sort by degree
        std::sort(last_level_s.begin(), last_level_s.end(), [this](const usize a, const usize b) {
            return degree[a] < degree[b];
        });

        /// Remove vertices that are *neighbours of vertices already considered* and get first 5 vertices only
        last_level.clear();
        for (const auto i: last_level_s) {
            bool discard = false;
            for (auto j_idx = row_index[i]; j_idx < row_index[i + 1]; ++j_idx) {
                const auto j = col_index[j_idx];

                if (considered_vertices.count(j) > 0) {
                    discard = true;
                    break;
                }
            }

            // Discard because the neighbor has already been visited
            if (discard) continue;

            last_level.push_back(i);
            considered_vertices.insert(i);
            
            // Limit to at most 5 elements
            if (last_level.size() >= 5)
                break; 
        }

        /// Set w(e) = infinity
        width_e = std::numeric_limits<usize>::max();

        /// (Test for termination) For each node i in f in order of ascending degree generate root level
        for (const usize i : last_level) {
            auto [eccentricity_i, width_i] = getEccentricityNWidth(i);

            // If h(i) > h(s) and w(i) < w(e) then set s = i and continue
            // Else, if w(i) < w(e) then set e = i and w(e) = w(i)
            if (eccentricity_i > eccentricity_s && width_i < width_e) {
                s = i;
                end = false;
                break; // continue
            } else if (width_i < width_e) {
                e = i;
                width_e = width_i;
            }

        }
    }

    /// Enhanced Sloan
    //Once we have a pseudoperipheral pair, we chose the node whose rooted level set has the greater depth or, if they have the same depth, the lesser width.
    auto [ecc_s, w_s] = getEccentricityNWidth(s);
    auto [ecc_e, w_e] = getEccentricityNWidth(e);
    if (ecc_e > ecc_s) {
        std::swap(s, e);
    } else if (ecc_e == ecc_s && w_e < w_s) {
        std::swap(s, e);
    }

    return {s, e};
}