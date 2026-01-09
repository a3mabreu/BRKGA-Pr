#include "csr.hpp"

CSR::CSR(usize rows, usize nnz) : m(rows), n_nz(nnz) {
    row_index.reserve(rows + 1);
    col_index.reserve(nnz);
}

CSR::CSR(const std::string& path, bool f_symmetric) {
    std::ifstream file(path);
    if (!file.is_open() || path.substr(path.length() - 4) != ".mtx") {
        std::cerr << "\n Failed to open file. Check the path!\n path: " << path << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    usize n_rows = 0, n_lines = 0;
    usize n_columns = 0, lines_read = 0;
    bool read_info = false;
    Element el;
    std::vector<Element> element_list;

    while (getline(file, line)) {
        if (n_lines == 0) {
            // Read information on first line
            if (!read_info) {
                std::istringstream iss(line);
                std::string last_word;
                // Read all words, keep the last one
                while (iss >> last_word) {};
                if (last_word.empty()) {
                    std::cerr << "Unable to read matriz info.";
                    exit(EXIT_FAILURE);
                }
                // Considering every matrix as symmetric
                if (last_word == "symmetric" || last_word == "skew-symmetric")
                    symmetric = true;
                else
                    symmetric = false;
                read_info = true;
                continue;
            }

            // Ignore comments
            if (line.empty() || line[0] == '%')
                continue;
            // Parse first line of file => (rows:m, columns:n, entries)
            std::istringstream iss(line);
            if (!(iss >> n_rows >> n_columns >> n_lines) || n_rows == 0 || n_lines == 0) {
                std::cerr << "Invalid header format\n";
                exit(EXIT_FAILURE);
            }
            if (n_rows != n_columns) {
                std::cerr << "m != n \n";
                exit(EXIT_FAILURE);
            }

            m = n_rows;
        } else {
            // Ignore comments
            if (line.empty() || line[0] == '%')
                continue;
            // Format => I1 J1 M(I1, J1)
            std::istringstream iss(line);
            if (!(iss >> el.i >> el.j)) {
                std::cerr << "Invalid line format\n";
                exit(EXIT_FAILURE);
            }
            std::string val;
            if ((iss >> val)) {
                el.v = stod(val);
            } else {
                // Default value when val is absent
                el.v = 0;
            }
            lines_read++;
            // Skip diagonal elements
            if (el.i == el.j) continue;

            // In these cases of *symmetric, skew-symmetric and Hermitian*, *only entries in the lower triangular portion need be supplied*
            if (symmetric && el.j > el.i) {
                std::cerr << "\nMTX Format Error: j > i. For symmetric matrix *only entries in the lower triangular portion need be supplied*\n";
                exit(EXIT_FAILURE);
            }

            // Adjust from 1-based to 0-based indexing
            el.i -= 1;
            el.j -= 1;

            element_list.push_back(el);
            // Add symmetric element if it's a symmetric matrix
            if (symmetric || f_symmetric) {
                Element el2;
                el2.i = el.j;
                el2.j = el.i;
                element_list.push_back(el2);
            }
        }
    }

    if (f_symmetric)
        symmetric = true;

    // Safety checks
    if (n_lines != lines_read) {
        std::cerr << "Unexpected number of lines read\n";
        exit(EXIT_FAILURE);
    }

    // Sort element_list by row and column
    sort(element_list.begin(), element_list.end(), [](const Element& a, const Element& b) {
        if (a.i < b.i) return true;
        if (a.i > b.i) return false;
        return a.j < b.j;
    });
    
    // Remove duplicates
    auto last = std::unique(element_list.begin(), element_list.end(), [](const Element& a, const Element& b) {
        return a.i == b.i && a.j == b.j;
    });
    element_list.erase(last, element_list.end());

    // Populate CSR row_index and col_index based on sorted elements
    col_index.resize(element_list.size());
    row_index.resize(m + 1);
    n_nz = element_list.size();  

    // Variables to track row indexing
    usize current_row = 0;
    row_index[0] = 0;

    for (usize i = 0; i < element_list.size(); i++) {
        const Element& e = element_list[i];
        // Column index
        col_index[i] = e.j;

        // If current element's row index is greater than the current row
        while (current_row < e.i) {
            row_index[current_row + 1] = i;
            current_row++;
        }
    }

    // Fill remaining row_index entries
    for (usize i = current_row + 1; i <= m; i++) {
        row_index[i] = n_nz;
    }

    // Fill Labels
    labels.resize(m);
    iota(labels.begin(), labels.end(), 0);

    /// Degrees
    degree.resize(m);
    min_degree = std::numeric_limits<usize>::max();
    for (usize i = 0; i < m; i++) {
        degree[i] = row_index[i + 1] - row_index[i];
        if (degree[i] > max_degree)
            max_degree = degree[i];

        if (degree[i] < min_degree)
            min_degree = degree[i];
    }

    distances.resize(m);
    visited.resize(m);
}


// Breadth-First Search (BFS) 
// Modified for finding the distances of all vertices to v
void CSR::bfs(usize v) {
    visited.assign(m, false);
    distances.assign(m, 0);
    visited[v] = true;
    usize dist = 0;
    std::queue<usize> q;

    q.push(v);
    while (!q.empty()) {
        const usize level_size = q.size();
        ++dist;
        for (usize i = 0; i < level_size; i++) {
            const auto u = q.front();
            q.pop();
            // Neighbors of u
            for (usize j = row_index[u]; j < row_index[u + 1]; ++j) {
                const usize w = col_index[j];
                if (!visited[w]) {
                    q.push(w);
                    visited[w] = true;
                    distances[w] = dist;
                }
            }
        }
    }
}

// Evaluate profile
void CSR::evaluateProfile() {
    usize small_neighbor_label;
    profile = 0;

    for (usize i = 0; i < m; i++) {
        const usize li = labels[i];
        if (li == 0) continue;  // if the LABEL is 0
        small_neighbor_label = li;
        /// For each neighbor of i
        for (usize j_idx = row_index[i]; j_idx < row_index[i + 1]; j_idx++) {
            const usize j = col_index[j_idx];
            const usize lj = labels[j];

            small_neighbor_label = std::min(small_neighbor_label, lj);
            // Break if it's the smallest possible
            if (small_neighbor_label == 0)
                break;
        }
        profile += li - small_neighbor_label;
    }

    if (profile < best_profile) 
        best_profile = profile;
}

// Get vertices from the last level structure and eccentricity
std::pair<std::vector<usize>, usize> CSR::getLastLevelAndEccentricity(usize v) {
    visited.assign(m, false);
    visited[v] = true;
    std::vector<usize> last_level;
    std::queue<usize> q;
    usize eccentricity = 0;

    q.push(v);
    while (!q.empty()) {
        const usize level_size = q.size();
        last_level.clear();
        ++eccentricity;
        for (usize i = 0; i < level_size; i++) {
            const auto u = q.front();
            q.pop();
            last_level.emplace_back(u);
            // Neighbors of u
            for (usize j = row_index[u]; j < row_index[u+1]; j++) {
                const usize w = col_index[j];
                if (!visited[w]) {
                    q.push(w);
                    visited[w] = true;
                }
            }
        }
    }
    --eccentricity; // Adjust for the +1 of first vertex
    return {std::move(last_level), eccentricity};
}

// Get eccentricity and level width (max) of v
std::pair<usize, usize> CSR::getEccentricityNWidth(usize v) {
    visited.assign(m, false);
    visited[v] = true;
    std::queue<usize> q;
    usize eccentricity = 0;
    usize width = 0;

    q.push(v);
    while (!q.empty()) {
        const usize level_size = q.size();
        width = std::max(width, level_size);
        ++eccentricity;
        for (usize i = 0; i < level_size; i++) {
            const auto u = q.front();
            q.pop();
            // Neighbors of u
            for (usize j = row_index[u]; j < row_index[u+1]; j++) {
                const usize w = col_index[j];
                if (!visited[w]) {
                    q.push(w);
                    visited[w] = true;
                }
            }
        }
    }
    --eccentricity; // Adjust for the +1 of first vertex

    return std::make_pair(eccentricity, width);
}

// Get the diameter of the graph (max eccentricity)
usize CSR::getDiameter() {
    // Eccentricity for all vertices
    std::vector<usize> eccentricity(m, 0);
    for (usize i = 0; i < m; i++) {
        // Find distances from i
        bfs(i);
        eccentricity[i] = *max_element(distances.begin(), distances.end());
    }

    // Returns maximum eccentricity
    return *max_element(eccentricity.begin(), eccentricity.end());
}

bool CSR::isFeasible() const {
    std::set<usize> unique_elements(labels.begin(), labels.end());
    
    if (unique_elements.size() != labels.size()) {
        std::cerr << "\n Labels are not unique!" << std::endl;
        return false;
    }        
    
    if (*unique_elements.begin() != 0) {
        std::cerr << "\n Labels do not start from 0!" << std::endl;
        return false;
    }

    if (*prev(unique_elements.end()) != (m - 1)) {
        std::cerr << "\n Labels do not end in m - 1!" << std::endl;
        return false;
    }
    
    usize l = 0;
    for (auto it = unique_elements.begin(); it != unique_elements.end(); ++it) {
        if (*it != l) {
            std::cerr << "\n Labels are not continuous!" << std::endl;
            return false;
        }
        ++l;
    }

    return true;
}
