#ifdef OCTOTIGER_HAVE_CUDA
#include "octotiger/multipole_interactions/calculate_stencil.hpp"
#include "octotiger/multipole_interactions/cuda_multipole_interaction_interface.hpp"
#include "octotiger/multipole_interactions/multipole_cuda_kernel.hpp"

#include "octotiger/defs.hpp"
#include "octotiger/options.hpp"

#include <array>
#include <vector>

namespace octotiger { namespace fmm { namespace multipole_interactions {

    cuda_multipole_interaction_interface::cuda_multipole_interaction_interface(
        void)
      : multipole_interaction_interface()
      , theta(opts().theta)
    {
    }

    void cuda_multipole_interaction_interface::compute_multipole_interactions(
        std::vector<real>& monopoles, std::vector<multipole>& M_ptr,
        std::vector<std::shared_ptr<std::vector<space_vector>>>& com_ptr,
        std::vector<neighbor_gravity_type>& neighbors, gsolve_type type,
        real dx, std::array<bool, geo::direction::count()>& is_direction_empty,
        std::array<real, NDIM> xbase)
    {
        kernel_scheduler::scheduler().init();
        // Check where we want to run this:
        int slot = kernel_scheduler::scheduler().get_launch_slot();
        if (slot == -1 || m2m_type == interaction_kernel_type::OLD)
        {
            // Run fallback CPU implementation
            multipole_interaction_interface::compute_multipole_interactions(
                monopoles, M_ptr, com_ptr, neighbors, type, dx,
                is_direction_empty, xbase);
        }
        else
        {    // run on cuda device
            if (type == RHO)
                cuda_launch_counter()++;
            else
                cuda_launch_counter_non_rho()++;
            // Move data into SoA arrays
            auto staging_area =
                kernel_scheduler::scheduler().get_staging_area(slot);
            update_input(monopoles, M_ptr, com_ptr, neighbors, type, dx, xbase,
                staging_area.local_monopoles, staging_area.local_expansions_SoA,
                staging_area.center_of_masses_SoA);

            // Queue moving of input data to device
            util::cuda_helper& gpu_interface =
                kernel_scheduler::scheduler().get_launch_interface(slot);
            kernel_device_enviroment& env =
                kernel_scheduler::scheduler().get_device_enviroment(slot);
            gpu_interface.copy_async(env.device_local_monopoles,
                staging_area.local_monopoles.data(), local_monopoles_size,
                cudaMemcpyHostToDevice);
            gpu_interface.copy_async(env.device_local_expansions,
                staging_area.local_expansions_SoA.get_pod(),
                local_expansions_size, cudaMemcpyHostToDevice);
            gpu_interface.copy_async(env.device_center_of_masses,
                staging_area.center_of_masses_SoA.get_pod(),
                center_of_masses_size, cudaMemcpyHostToDevice);

            // Launch kernel and queue copying of results
            dim3 const grid_spec(INX, 1, 1);
            dim3 const threads_per_block(1, INX, INX);
            if (type == RHO)
            {
                bool second_phase = false;
                void* args[] = {&(env.device_local_monopoles),
                    &(env.device_center_of_masses),
                    &(env.device_local_expansions),
                    &(env.device_potential_expansions),
                    &(env.device_angular_corrections), &theta, &second_phase};
                gpu_interface.execute(
                    (void const*) &cuda_multipole_interactions_kernel_rho,
                    grid_spec, threads_per_block, args, 0);
                // second_phase = true;
                // gpu_interface.execute(&cuda_multipole_interactions_kernel_rho, grid_spec,
                //                       threads_per_block, args, 0);
                gpu_interface.copy_async(angular_corrections_SoA.get_pod(),
                    env.device_angular_corrections, angular_corrections_size,
                    cudaMemcpyDeviceToHost);
            }
            else
            {
                bool second_phase = false;
                void* args[] = {&(env.device_local_monopoles),
                    &(env.device_center_of_masses),
                    &(env.device_local_expansions),
                    &(env.device_potential_expansions), &theta, &second_phase};
                gpu_interface.execute(
                    (void const*) &cuda_multipole_interactions_kernel_non_rho,
                    grid_spec, threads_per_block, args, 0);
                // second_phase = true;
                // gpu_interface.execute(&cuda_multipole_interactions_kernel_non_rho, grid_spec,
                //     threads_per_block, args, 0);
            }
            gpu_interface.copy_async(potential_expansions_SoA.get_pod(),
                env.device_potential_expansions, potential_expansions_size,
                cudaMemcpyDeviceToHost);

            // Wait for stream to finish and allow thread to jump away in the meantime
            auto fut = gpu_interface.get_future();
            fut.get();

            // Copy results back into non-SoA array
            potential_expansions_SoA.add_to_non_SoA(grid_ptr->get_L());
            if (type == RHO)
                angular_corrections_SoA.to_non_SoA(grid_ptr->get_L_c());
        }
    }

}}}
#endif
