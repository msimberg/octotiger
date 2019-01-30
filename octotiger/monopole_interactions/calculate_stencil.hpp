#pragma once

#include "octotiger/common_kernel/multiindex.hpp"

#include "octotiger/real.hpp"

#include <array>
#include <utility>
#include <octotiger/debug_vector.hpp>

namespace octotiger {
namespace fmm {
    namespace monopole_interactions {

        std::pair<oct::vector<multiindex<>>, oct::vector<oct::array<real, 4>>> calculate_stencil();
        std::pair<oct::vector<bool>, oct::vector<oct::array<real, 4>>>
        calculate_stencil_masks(oct::vector<multiindex<>> superimposed_stencil);

    }    // namespace monopole_interactions
}    // namespace fmm
}    // namespace octotiger
