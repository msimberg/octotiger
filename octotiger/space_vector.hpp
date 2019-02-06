/*
 * space_vector_gen.hpp
 *
 *  Created on: Jun 3, 2015
 *      Author: dmarce1
 */

#ifndef SPACE_VECTOR_HPP_
#define SPACE_VECTOR_HPP_

#include "octotiger/defs.hpp"
#include "octotiger/real.hpp"

#include <hpx/parallel/traits/vector_pack_type.hpp>
#include <hpx/runtime/serialization/datapar.hpp>
#include <hpx/traits/is_bitwise_serializable.hpp>

#include <Vc/Vc>

// #if !defined(HPX_HAVE_DATAPAR)


// #else

// #if defined(Vc_HAVE_AVX512F) || defined(Vc_HAVE_AVX)

#ifdef __AVX2__
using space_vector = Vc::Vector<real, Vc::VectorAbi::Avx>;
#else
using space_vector = hpx::parallel::traits::vector_pack_type<real, 4>::type;
#endif

#endif /* SPACE_VECTOR_HPP_ */
