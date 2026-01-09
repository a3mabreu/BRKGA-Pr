#ifndef MAX_HEAP_RH_H
#define MAX_HEAP_RH_H

#include "types.hpp"

//// Everething in inlined

// Get first (max element), do max‑heapify (down) and update idx_H
inline VertexCost getFirst(std::vector<VertexCost>& H, robin_hood::unordered_map<usize, usize>& idx_H) {
    const auto u = H.front();

    // Remove the first element from idx_H
    idx_H.erase(u.first);

    // Replace the root with the last element in the heap
    H[0] = H.back();
    H.pop_back();

    if (H.empty()) {
        return u; // Heap is now empty
    }

    // Update idx_H for the new root
    idx_H[H[0].first] = 0;

    // Max-heapify from the root (bubble down)
    usize heap_idx = 0;
    const usize last_idx = H.size() - 1;
    while (true) {
        const usize left_child = 2 * heap_idx + 1;
        const usize right_child = 2 * heap_idx + 2;
        usize largest = heap_idx;

        // Find the largest child
        if (left_child <= last_idx && (H[left_child].second > H[largest].second)) {
            largest = left_child;
        }
        if (right_child <= last_idx && (H[right_child].second > H[largest].second)) {
            largest = right_child;
        }

        if (largest != heap_idx) {
            // Swap with the largest child
            std::swap(H[heap_idx], H[largest]);
            // Update idx_H
            idx_H[H[heap_idx].first] = heap_idx;
            idx_H[H[largest].first] = largest;
            // Move down to the swapped position
            heap_idx = largest;
        } else {
            break; // Heap property is satisfied
        }
    }

    return u;
}

// Max-heapify UP (bubble up) and update idx_H
// This function is used when a key is increased.
inline void bubleUp(std::vector<VertexCost>& H, robin_hood::unordered_map<usize, usize>& idx_H, usize w, int weight) {
    // Update cost in max heap
    H[idx_H.at(w)].second = weight;
    usize heap_idx = idx_H.at(w);

    // Restore heap property by moving up: parent's value must be >= child's value.
    while (heap_idx > 0) {
        const usize parent_idx = (heap_idx - 1) / 2;
        // If parent's value is already larger, we're done.
        if (H[parent_idx].second >= H[heap_idx].second) {
            break; 
        }
        
        // Swap parent and current node in the heap
        std::swap(H[parent_idx], H[heap_idx]);

        // Update idx_H
        idx_H[H[parent_idx].first] = parent_idx;
        idx_H[H[heap_idx].first] = heap_idx;

        // Move up to the parent node
        heap_idx = parent_idx;
    }
}


// Max-heapify DOWN (bubble down) and update idx_H
// This function is used when a key is decreased.
inline void bubbleDown(std::vector<VertexCost>& H, robin_hood::unordered_map<usize, usize>& idx_H, usize w, int weight) {
    // Update cost in max heap
    usize heap_idx = idx_H.at(w);
    H[heap_idx].second = weight;
    const usize size = H.size();

    // Restore heap property by moving down:
    // Each parent must be greater than or equal to its children.
    while (true) {
        const usize left_child = 2 * heap_idx + 1;
        const usize right_child = 2 * heap_idx + 2;
        usize largest = heap_idx;

        // Find the largest among the current node and its children
        if (left_child < size && H[left_child].second > H[largest].second) {
            largest = left_child;
        }
        if (right_child < size && H[right_child].second > H[largest].second) {
            largest = right_child;
        }

        // If the largest is still the current node, stop
        if (largest == heap_idx) {
            break;
        }

        // Swap the current node with the largest child
        std::swap(H[heap_idx], H[largest]);

        // Update idx_H
        idx_H[H[heap_idx].first] = heap_idx;
        idx_H[H[largest].first] = largest;

        // Move down to the swapped position
        heap_idx = largest;
    }
}


inline void insertHeap(std::vector<VertexCost>& H, robin_hood::unordered_map<usize, usize>& idx_H, usize vertex, int cost) {
    const usize size = H.size();
    // Add new element at the end of the heap
    H.emplace_back(vertex, cost);
    // Store its index in idx_H
    idx_H[vertex] = size;
    // Bubble up to maintain heap property
    bubleUp(H, idx_H, vertex, cost);
}


// Remove an element from the heap H and idx_H and adjust the heap
inline void removeElement(std::vector<VertexCost>& H, robin_hood::unordered_map<usize, usize>& idx_H, usize w) {
    // Check if the element exists in the heap
    // if (idx_H.find(w) == idx_H.end()) return;

    const usize heap_idx = idx_H.at(w);  // Get element index
    const usize last_idx = H.size() - 1;

    // If it's the last element, simply pop it
    if (heap_idx == last_idx) {
        idx_H.erase(w);
        H.pop_back();
        return;
    }

    // Replace the removed element with the last element
    H[heap_idx] = H[last_idx];
    idx_H[H[heap_idx].first] = heap_idx;  // Update idx_H
    H.pop_back();
    idx_H.erase(w);  // Remove w from idx_H

    // Restore heap property:
    // If the new value is larger than its parent, bubble up;
    // otherwise, bubble down.
    if (heap_idx > 0 && H[heap_idx].second > H[(heap_idx - 1) / 2].second) {
        bubleUp(H, idx_H, H[heap_idx].first, H[heap_idx].second);
    } else {
        bubbleDown(H, idx_H, H[heap_idx].first, H[heap_idx].second);
    }
}


#endif
