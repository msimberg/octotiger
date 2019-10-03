#if !defined(__CUDA_ARCH__)
#include "kernel.hpp"

#include <hpx/hpx_main.hpp>
#include <Kokkos_Core.hpp>

#include <iostream>

int main(int argc, char** argv) {
    printf("hello from hpx_main\n");

    Kokkos::initialize(argc, argv);
    Kokkos::print_configuration(std::cout);
    kernel();
    Kokkos::finalize();

    return 0;
}
#endif
