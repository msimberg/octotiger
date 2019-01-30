#include "octotiger/multipole_interactions/multipole_interaction_interface.hpp"
#include "octotiger/multipole_interactions/calculate_stencil.hpp"
#include "octotiger/multipole_interactions/multipole_cpu_kernel.hpp"

#include "octotiger/common_kernel/interactions_iterators.hpp"

#include "octotiger/options.hpp"

#include <algorithm>
#include <array>
#include <octotiger/debug_vector.hpp>

// Big picture questions:
// - use any kind of tiling?

namespace octotiger {
namespace fmm {
    namespace multipole_interactions {

        thread_local const two_phase_stencil multipole_interaction_interface::stencil =
            calculate_stencil();
        thread_local oct::vector<real>
            multipole_interaction_interface::local_monopoles_staging_area(EXPANSION_COUNT_PADDED);
        thread_local struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>
        multipole_interaction_interface::local_expansions_staging_area;
        thread_local struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>
        multipole_interaction_interface::center_of_masses_staging_area;
        thread_local const oct::vector<bool> multipole_interaction_interface::stencil_masks =
            calculate_stencil_masks(multipole_interaction_interface::stencil).first;
        thread_local const oct::vector<bool> multipole_interaction_interface::inner_stencil_masks =
            calculate_stencil_masks(multipole_interaction_interface::stencil).second;


        multipole_interaction_interface::multipole_interaction_interface(void) {
            local_monopoles_staging_area = oct::vector<real>(ENTRIES);
            this->m2m_type = opts().m2m_kernel_type;
        }

        void multipole_interaction_interface::compute_multipole_interactions(
            oct::vector<real>& monopoles, oct::vector<multipole>& M_ptr,
            oct::vector<std::shared_ptr<oct::vector<space_vector>>>& com_ptr,
            oct::vector<neighbor_gravity_type>& neighbors, gsolve_type type, real dx,
            oct::array<bool, geo::direction::count()>& is_direction_empty,
            oct::array<real, NDIM> xbase) {
            update_input(monopoles, M_ptr, com_ptr, neighbors, type, dx, xbase,
                local_monopoles_staging_area, local_expansions_staging_area,
                center_of_masses_staging_area);
            compute_interactions(is_direction_empty, neighbors, local_monopoles_staging_area,
                local_expansions_staging_area, center_of_masses_staging_area);
        }

        void multipole_interaction_interface::compute_interactions(
            oct::array<bool, geo::direction::count()>& is_direction_empty,
            oct::vector<neighbor_gravity_type>& all_neighbor_interaction_data,
            const oct::vector<real>& local_monopoles,
            const struct_of_array_data<expansion, real, 20, ENTRIES, SOA_PADDING>&
                local_expansions_SoA,
            const struct_of_array_data<space_vector, real, 3, ENTRIES, SOA_PADDING>&
                center_of_masses_SoA) {
            if (m2m_type == interaction_kernel_type::SOA_CPU) {
                struct_of_array_data<expansion, real, 20, INNER_CELLS, SOA_PADDING>
                    potential_expansions_SoA;
                struct_of_array_data<space_vector, real, 3, INNER_CELLS, SOA_PADDING>
                    angular_corrections_SoA;

                multipole_cpu_kernel kernel;
                kernel.apply_stencil_non_blocked(local_expansions_SoA, center_of_masses_SoA,
                    potential_expansions_SoA, angular_corrections_SoA, local_monopoles,
                                                 stencil_masks, inner_stencil_masks, type);

                if (type == RHO) {
                    angular_corrections_SoA.to_non_SoA(grid_ptr->get_L_c());
                }
                potential_expansions_SoA.add_to_non_SoA(grid_ptr->get_L());
            } else {
                // old-style interaction calculation
                // computes inner interactions
                grid_ptr->compute_interactions(type);
                // computes boundary interactions
                for (auto const& dir : geo::direction::full_set()) {
                    if (!is_direction_empty[dir]) {
                        neighbor_gravity_type& neighbor_data = all_neighbor_interaction_data[dir];
                        grid_ptr->compute_boundary_interactions(type, neighbor_data.direction,
                            neighbor_data.is_monopole, neighbor_data.data);
                    }
                }
            }
        }
    }    // namespace multipole_interactions
}    // namespace fmm
}    // namespace octotiger
