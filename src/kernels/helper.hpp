#pragma once

#include "simd.hpp"

namespace octotiger {
namespace fmm {
    namespace detail {

        // calculates 1/distance between i and j
        inline real reciprocal_distance(const integer i0, const integer i1, const integer i2,
            const integer j0, const integer j1, const integer j2) {
            real tmp = (sqr(i0 - j0) + sqr(i1 - j1) + sqr(i2 - j2));
            // protect against sqrt(0)
            if (tmp > 0.0) {                      // TODO: remove this branch
                return 1.0 / (std::sqrt(tmp));    // TODO: use squared theta instead for comparison
            } else {
                return 1.0e+10;    // dummy value
            }
        }

        // calculates 1/distance between i and j
        // TODO: change to simd_vector
        inline real distance_squared(const multiindex<>& i, const multiindex<>& j) {
            // protect against sqrt(0)
            // return 1.0 / (std::sqrt(tmp));    // TODO: use squared theta instead for compariso
            return (sqr(i.x - j.x) + sqr(i.y - j.y) + sqr(i.z - j.z));
        }

        // calculates 1/distance between i and j
        // TODO: change to simd_vector
        inline int_simd_vector distance_squared(
            const multiindex<int_simd_vector>& i, const multiindex<int_simd_vector>& j) {
            // protect against sqrt(0)
            // return 1.0 / (std::sqrt(tmp));
            // TODO: use squared theta instead for compariso
            return (sqr(i.x - j.x) + sqr(i.y - j.y) + sqr(i.z - j.z));
        }

        // checks whether the index tuple is inside the current node
        inline bool is_interior(const integer i0, const integer i1, const integer i2) {
            bool check = true;
            if (i0 < 0 || i0 >= G_NX) {
                check = false;
            } else if (i1 < 0 || i1 >= G_NX) {
                check = false;
            } else if (i2 < 0 || i2 >= G_NX) {
                check = false;
            }
            return check;
        }
    }    // namespace detail
}    // namespace fmm
}    // namespace octotiger