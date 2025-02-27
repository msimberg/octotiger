#pragma once

#include "octotiger/defs.hpp"

#include <cinttypes>
#include <cstddef>
#include <cstdint>

namespace octotiger {

namespace detail {

    constexpr uint64_t const_pow(uint64_t base, uint64_t exp) {
        return exp >= 1 ? base * const_pow(base, exp - 1) : 1;
    }
}

namespace fmm {

    // constants from defs.hpp in deobfuscated form
    constexpr uint64_t DIMENSION = NDIM;    // 3

    // number of expansions in each cell per direction
    constexpr int64_t INNER_CELLS_PER_DIRECTION = INX;    // 8

    // number of expansions in each cell (inner cells)
    constexpr uint64_t INNER_CELLS = detail::const_pow(INNER_CELLS_PER_DIRECTION, DIMENSION);

    constexpr uint64_t INNER_CELLS_PADDING_DEPTH = INX;
    constexpr uint64_t PADDING_OFFSET = INNER_CELLS_PADDING_DEPTH - 5;

constexpr uint64_t PADDED_STRIDE = INNER_CELLS_PER_DIRECTION + 2 * (5);

    constexpr uint64_t ENTRIES = PADDED_STRIDE * PADDED_STRIDE * PADDED_STRIDE;

    constexpr uint64_t EXPANSION_COUNT_PADDED = detail::const_pow(PADDED_STRIDE, DIMENSION);
    constexpr uint64_t EXPANSION_COUNT_NOT_PADDED = INNER_CELLS;

    constexpr uint64_t SOA_PADDING = 19;    // to prevent some of the 4k aliasing

    constexpr uint64_t STENCIL_SIZE = 1074;
    constexpr uint64_t STENCIL_BLOCKING = 16;
    constexpr int STENCIL_INX = 11;
    constexpr int STENCIL_MIN = -5;
    constexpr int STENCIL_MAX = 5;
    constexpr int FULL_STENCIL_SIZE = STENCIL_INX * STENCIL_INX * STENCIL_INX;

//constexpr uint64_t STENCIL_SIZE = 982;
    constexpr size_t NUMBER_LOCAL_MONOPOLE_VALUES = 1 * (ENTRIES);
    constexpr size_t NUMBER_LOCAL_EXPANSION_VALUES = 20 * (ENTRIES + SOA_PADDING);
    constexpr size_t NUMBER_MASS_VALUES = 3 * (ENTRIES + SOA_PADDING);
    constexpr size_t NUMBER_POT_EXPANSIONS = 20 * (INNER_CELLS + SOA_PADDING);
    constexpr size_t NUMBER_POT_EXPANSIONS_SMALL = 4 * (INNER_CELLS + SOA_PADDING);
    constexpr size_t NUMBER_ANG_CORRECTIONS = 3 * (INNER_CELLS + SOA_PADDING);
    constexpr size_t NUMBER_FACTORS = 20;

    constexpr uint64_t P2P_CHUNKSIZE = 128;
    constexpr uint64_t P2P_CHUNK_STENCIL_SIZE = P2P_CHUNKSIZE + 1;
    constexpr uint64_t P2P_PADDED_STENCIL_SIZE = STENCIL_SIZE + P2P_CHUNK_STENCIL_SIZE - STENCIL_SIZE % P2P_CHUNK_STENCIL_SIZE;

}    // namespace fmm
}    // namespace octotiger
