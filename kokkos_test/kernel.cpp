#include "kernel.hpp"

#include <Kokkos_Core.hpp>

#include <iostream>

void kernel() {
    int const n = 10;
    Kokkos::View<double*> a("a", n);
    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n), KOKKOS_LAMBDA(int const i) { printf("hello from default space\n"); a(i) = i; });
    Kokkos::fence();
    printf("DefaultExecutionSpace done\n");
    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n), KOKKOS_LAMBDA(int const i) { printf("hello from default host space\n"); });
    Kokkos::fence();
    printf("DefaultHostExecutionSpace done\n");
}
