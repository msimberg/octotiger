/*
 * rad_grid.hpp
 *
 *  Created on: May 20, 2016
 *      Author: dmarce1
 */

#ifndef RAD_GRID_HPP_
#define RAD_GRID_HPP_

#include "octotiger/defs.hpp"
#include "octotiger/geometry.hpp"
#include "octotiger/physcon.hpp"
#include "octotiger/real.hpp"
#include "octotiger/silo.hpp"
//#include "octotiger/sphere_points.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <octotiger/debug_vector.hpp>

#define R_BW 3
#define R_NX (INX+2*R_BW)
#define R_N3 (R_NX*R_NX*R_NX)

//using real = real;

class rad_grid {
public:
	static constexpr integer er_i = 0;
	static constexpr integer fx_i = 1;
	static constexpr integer fy_i = 2;
	static constexpr integer fz_i = 3;
private:
	static constexpr integer DX = R_NX * R_NX;
	static constexpr integer DY = R_NX;
	static constexpr integer DZ = 1;
	static std::unordered_map<std::string, int> str_to_index;
	static std::unordered_map<int, std::string> index_to_str;
	real dx;
	oct::array<oct::vector<real>, NRF> U;
	oct::array<oct::vector<real>, NRF> U0;
	oct::array<oct::array<oct::vector<real>, NRF>,NDIM> flux;
	oct::array<oct::array<oct::vector<real>*, NDIM>, NDIM> P;
	oct::vector<oct::vector<real>> X;
	oct::vector<real> mmw, X_spc, Z_spc;
	void reconstruct(oct::array<oct::vector<real>, NRF>&,oct::array<oct::vector<real>, NRF>&,int dir);
public:
	static void static_init();
	static oct::vector<std::string> get_field_names();
	void set(const std::string name, real* data);
	oct::vector<silo_var_t> var_data() const;
	void set_X( const oct::vector<oct::vector<real>>& x );
	void restore();
	void store();

	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & dx;
		arc & U;
	}
	void compute_mmw(const oct::vector<oct::vector<real>>& U);
	void change_units(real m, real l, real t, real k);
	real rad_imp_comoving(real& E, real& e, real rho, real mmw, real X, real Z, real dt);
	void sanity_check();
	void initialize_erad(const oct::vector<real> rho, const oct::vector<real> tau);
	void set_dx(real dx);
	//void compute_fEdd();
	void compute_flux();
	void advance(real dt, real beta);
	void rad_imp(oct::vector<real>& egas, oct::vector<real>& tau, oct::vector<real>& sx, oct::vector<real>& sy, oct::vector<real>& sz,
			const oct::vector<real>& rho, real dt);
	oct::vector<real> get_restrict() const;
	oct::vector<real> get_prolong(const oct::array<integer, NDIM>& lb, const oct::array<integer, NDIM>& ub);
	void set_prolong(const oct::vector<real>&);
	void set_restrict(const oct::vector<real>&, const geo::octant&);
	void set_flux_restrict(const oct::vector<real>& data, const oct::array<integer, NDIM>& lb, const oct::array<integer, NDIM>& ub,
			const geo::dimension& dim);
	oct::vector<real> get_flux_restrict(const oct::array<integer, NDIM>& lb, const oct::array<integer, NDIM>& ub, const geo::dimension& dim) const;
	oct::vector<real> get_intensity(const oct::array<integer, NDIM>& lb, const oct::array<integer, NDIM>& ub, const geo::octant&);
	void allocate();
	rad_grid(real dx);
	rad_grid();
	void set_boundary(const oct::vector<real>& data, const geo::direction& dir);
	real get_field(integer f, integer i, integer j, integer k) const;
	void set_field(real v, integer f, integer i, integer j, integer k);
	void set_physical_boundaries(geo::face f, real t);
	oct::vector<real> get_boundary(const geo::direction& dir);
	using kappa_type = std::function<real(real)>;

	real hydro_signal_speed(const oct::vector<real>& egas, const oct::vector<real>& tau, const oct::vector<real>& sx, const oct::vector<real>& sy, const oct::vector<real>& sz,
			const oct::vector<real>& rho);

	friend class node_server;
};




#endif /* RAD_GRID_HPP_ */

