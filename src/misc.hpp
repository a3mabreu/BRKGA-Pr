#ifndef MISC_H
#define MISC_H

#include "types.hpp"
#include <functional>
#include <random>

#ifdef SEED
  inline std::mt19937 mt(static_cast<std::mt19937::result_type>(SEED));
#else
  inline std::random_device rd;
  inline std::mt19937 mt(rd());
#endif
// Static distributions
constexpr realT after_one = std::nextafter(1.0, 2.0);
inline std::uniform_real_distribution<realT> disRK(0.0, 1.0);
inline std::uniform_int_distribution<usize> distRNumber;
inline std::uniform_int_distribution<usize> distRLabel;
inline std::uniform_real_distribution<realT> betweenZeroAndOne(0.0, after_one);

////////        INLINE
// Get Mersenne Twister generator
inline std::mt19937& getMT() {
    return mt;
}

// Generating an USIZE random number between min and max (inclusive)
inline usize usizeRandomNumber(const usize min, const usize max) {
    distRNumber.param(typename decltype(distRNumber)::param_type(min, max));
    return distRNumber(getMT());
}

// Generating an USIZE random number between [0, max] (inclusive)
// Do not work if multiple graphs are instantiated
inline usize usizeRandomLabel() {
    return distRLabel(getMT());
}

// Generating an REAL random number between [0, 1)
inline realT realRK() {
    return disRK(getMT());
}

// Generating an REAL random number between [0, 1] *Inclusive*
inline realT realZeroOneInclusive() {
    return betweenZeroAndOne(getMT());
}

#endif