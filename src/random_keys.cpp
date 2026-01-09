#include "random_keys.hpp"
#include "misc.hpp"

// Decode random keys into a solution to labels (SORT)
void decoder(CSR& csr) {
    const usize n = csr.m;
    
    for (usize i = 0; i < n; ++i)
        csr.indexed_rk[i] = {csr.random_keys[i], i};
    
    // Sort by the random key values
    std::sort(csr.indexed_rk.begin(), csr.indexed_rk.end(), [](const IndexedElement& a, const IndexedElement& b) {
        return a.value < b.value;
    });
    
    // Assign labels based on sorted order
    for (usize i = 0; i < n; ++i)
        csr.labels[i] = csr.indexed_rk[i].index;
}

// Encode a solution (Labels) into RK representation
// Warning !!! "realT" must be double (in types.hpp) for encoder to work
void encoder(CSR& csr) {
    const usize n = csr.m;

    // Generate RKs
    for (usize i = 0; i < n; ++i)
        csr.tmp_rk[i] = realRK();

    // Sort the random key values by labels
    std::sort(csr.tmp_rk.begin(), csr.tmp_rk.end());
    for (usize i = 0; i < n; ++i)
        csr.random_keys[csr.labels[i]] = csr.tmp_rk[i];
}