/*
 * silo.cpp
 *
 *  Created on: Oct 16, 2018
 *      Author: dmarce1
 */

#include "octotiger/silo.hpp"
#include "octotiger/node_registry.hpp"
#include "octotiger/node_server.hpp"
#include "octotiger/options.hpp"
#include "octotiger/physcon.hpp"
#include "octotiger/util.hpp"

#include <hpx/lcos/broadcast.hpp>
#include <hpx/util/format.hpp>
#include <hpx/util/io_service_pool.hpp>

#include <future>
#include <mutex>
#include <set>
#include <vector>

#define OUTPUT_ROCHE
template<class T>
struct read_silo_var {
	T operator()(DBfile* db, const char* name) const {
		int one = 1;
		T var;
		if (DBReadVar(db, name, &var) != 0) {
			std::cout << "Unable to read " << name << " \n";
		}
		return var;
	}
};


struct node_entry_t {
	bool load;
	integer position;
	integer locality_id;
	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & load;
		arc & position;
		arc & locality_id;
	}
};

using dir_map_type = std::unordered_map<node_location::node_id, node_entry_t>;
dir_map_type node_dir_;
DBfile* db_;

static const auto& localities = options::all_localities;

static std::mutex silo_mtx_;


#define SILO_TEST(i) \
	if( i != 0 ) printf( "SILO call failed at %i\n", __LINE__ );

void load_options_from_silo(std::string fname, DBfile* db) {
	const auto func =
			[&fname,&db]()
			{
				bool leaveopen;
				if( db == NULL )
				{
					db = DBOpenReal( fname.c_str(), SILO_DRIVER, DB_READ);
					leaveopen = false;
				}
				else
				{
					leaveopen = true;
				}
				if (db != NULL)
				{

					read_silo_var<integer> ri;
					read_silo_var<real> rr;
					integer version = ri( db, "version");
					if( version > SILO_VERSION) {
						printf( "WARNING: Attempting to load a version %i SILO file, maximum version allowed for this Octo-tiger is %i\n", int(version), SILO_VERSION);
					}
					opts().code_to_g = rr(db, "code_to_g");
					opts().code_to_s = rr(db, "code_to_s");
					opts().code_to_cm = rr(db, "code_to_cm");
					opts().n_species = ri(db, "n_species");
					opts().eos = eos_type(ri(db, "eos"));
					opts().gravity = ri(db, "gravity");
					opts().hydro = ri(db, "hydro");
					opts().omega = rr(db, "omega") * opts().code_to_s;
					opts().output_dt = rr(db, "output_frequency");
					opts().problem = problem_type(ri(db, "problem"));
					opts().radiation = ri(db, "radiation");
					opts().refinement_floor = rr(db, "refinement_floor");
					opts().xscale = rr(db, "xscale");
					opts().atomic_number.resize(opts().n_species);
					opts().atomic_mass.resize(opts().n_species);
					opts().X.resize(opts().n_species);
					opts().Z.resize(opts().n_species);
					SILO_TEST(DBReadVar(db, "atomic_number", opts().atomic_number.data()));
					SILO_TEST(DBReadVar(db, "atomic_mass", opts().atomic_mass.data()));
					SILO_TEST(DBReadVar(db, "X", opts().X.data()));
					SILO_TEST(DBReadVar(db, "Z", opts().Z.data()));
					if( !leaveopen)
					{
						DBClose(db);
					}
				}
				else
				{
					std::cout << "Could not load " << fname;
					throw;
				}
			};
	if (db == NULL) {
		GET(hpx::threads::run_as_os_thread(func));
	} else {
		func();
	}
	grid::set_omega(opts().omega, false);
	set_units(1./opts().code_to_g, 1./opts().code_to_cm, 1./opts().code_to_s, 1); /**/

}

extern real output_time;
extern real output_rotation_count;

void load_open(std::string fname, dir_map_type map) {
	printf("LOAD OPENED on proc %i\n", hpx::get_locality_id());
	load_options_from_silo(fname, db_); /**/
	hpx::threads::run_as_os_thread([&]() {
		db_ = DBOpenReal( fname.c_str(), SILO_DRIVER, DB_READ);
		read_silo_var<real> rr;
		output_time = rr(db_, "cgs_time"); /**/
		output_rotation_count = 2 * M_PI * rr(db_, "rotational_time"); /**/
		printf( "rotational_time = %e\n", output_rotation_count);
		output_time /= opts().code_to_s;
		node_dir_ = std::move(map);
		printf( "%e\n", output_time );
//		sleep(100);
		}).get();
}

void load_close() {
	DBClose(db_);
}

HPX_PLAIN_ACTION(load_close, load_close_action);
HPX_PLAIN_ACTION(load_open, load_open_action);

node_server::node_server(const node_location& loc) :
		my_location(loc) {
	const auto& localities = opts().all_localities;
	initialize(0.0, 0.0);
	step_num = gcycle = hcycle = rcycle = 0;

	auto iter = node_dir_.find(loc.to_id());
	assert(iter != node_dir_.end());

	if (!iter->second.load) {
		printf("Creating %s on %i\n", loc.to_str().c_str(), int(hpx::get_locality_id()));
		int nc = 0;
		for (int ci = 0; ci < NCHILD; ci++) {
			auto cloc = loc.get_child(ci);
			auto iter = node_dir_.find(cloc.to_id());
			if (iter != node_dir_.end() && iter->second.locality_id != hpx::get_locality_id()) {
				is_refined = true;
				children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc);
				nc++;
			}
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			auto cloc = loc.get_child(ci);
			auto iter = node_dir_.find(cloc.to_id());
			if (iter != node_dir_.end() && iter->second.locality_id == hpx::get_locality_id()) {
				is_refined = true;
				children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc);
				nc++;
			}
		}
		assert(nc == 0 || nc == NCHILD);
	} else {
		printf("Loading %s on %i\n", loc.to_str().c_str(), int(hpx::get_locality_id()));
		silo_load_t load;
		static const auto hydro_names = grid::get_hydro_field_names();
		load.vars.resize(hydro_names.size());
		load.outflows.resize(hydro_names.size());
		hpx::threads::run_as_os_thread([&]() {
			static std::mutex mtx;
			std::lock_guard<std::mutex> lock(mtx);
			const std::string suffix = oct_to_str(loc.to_id());
			for( int f = 0; f != hydro_names.size(); f++) {
				const auto this_name = std::string( "/") + suffix + std::string( "/") + hydro_names[f]; /**/
				auto var = DBGetQuadvar(db_,this_name.c_str());
				load.nx = var->dims[0];
				const int nvar = load.nx * load.nx * load.nx;
				load.outflows[f].first = load.vars[f].first = hydro_names[f];
				load.vars[f].second.resize(nvar);
				read_silo_var<real> rd;
				load.outflows[f].second = rd(db_, outflow_name(this_name).c_str());
				std::memcpy(load.vars[f].second.data(), var->vals[0], sizeof(real)*nvar);
				DBFreeQuadvar(var);
			}
		}).get();
		if (load.nx == INX) {
			is_refined = false;
			for (integer f = 0; f < hydro_names.size(); f++) {
				grid_ptr->set(load.vars[f].first, load.vars[f].second.data());
				grid_ptr->set_outflow(std::move(load.outflows[f]));
			}
			grid_ptr->rho_from_species();
		} else {
			is_refined = true;
			auto child_loads = load.decompress();
			for (integer ci = 0; ci < NCHILD; ci++) {
				auto cloc = loc.get_child(ci);
				auto iter = node_dir_.find(cloc.to_id());
				children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc, child_loads[ci]);
			}
		}
	}
	current_time = output_time;
	rotational_time = output_rotation_count;
}

node_server::node_server(const node_location& loc, silo_load_t load) :
		my_location(loc) {
	printf("Distributing %s on %i\n", loc.to_str().c_str(), int(hpx::get_locality_id()));
	const auto& localities = opts().all_localities;
	initialize(0.0, 0.0);
	step_num = gcycle = hcycle = rcycle = 0;
	int nc = 0;
	static const auto hydro_names = grid::get_hydro_field_names();
	if (load.nx == INX) {
		is_refined = false;
		for (integer f = 0; f < hydro_names.size(); f++) {
			grid_ptr->set(load.vars[f].first, load.vars[f].second.data());
			grid_ptr->set_outflow(std::move(load.outflows[f]));
		}
		grid_ptr->rho_from_species();
	} else {
		is_refined = true;
		auto child_loads = load.decompress();
		for (integer ci = 0; ci < NCHILD; ci++) {
			auto cloc = loc.get_child(ci);
			auto iter = node_dir_.find(cloc.to_id());
			children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc, child_loads[ci]);
		}
	}
	current_time = output_time;
	rotational_time = output_rotation_count;
	assert(nc == 0 || nc == NCHILD);
}

extern int epoch;

void load_data_from_silo(std::string fname, node_server* root_ptr, hpx::id_type root) {
	const integer nprocs = opts().all_localities.size();
	static int sz = localities.size();
	DBfile* db = GET(hpx::threads::run_as_os_thread(DBOpenReal, fname.c_str(), SILO_DRIVER, DB_READ));
	epoch = GET(hpx::threads::run_as_os_thread(read_silo_var<integer>(), db, "epoch"));
	epoch++;
	std::vector<node_location::node_id> node_list;
	std::vector<integer> positions;
	std::vector<hpx::future<void>> futs;
	int node_count;
	if (db != NULL) {
		DBmultimesh* master_mesh = GET(hpx::threads::run_as_os_thread([&]()
		{
			return DBGetMultimesh( db, "quadmesh");
		}));
		const int chunk_size = std::ceil(real(master_mesh->nblocks) / real(sz));
		hpx::threads::run_as_os_thread([&]() {
			const read_silo_var<integer> ri;
			node_count = ri(db,"node_count");
			node_list.resize(node_count);
			positions.resize(node_count);
			DBReadVar(db, "node_list", node_list.data());
			DBReadVar(db, "node_positions", positions.data());
		}).get();
		GET(hpx::threads::run_as_os_thread(DBClose, db));
		std::set<node_location::node_id> load_locs;
		for (int i = 0; i < master_mesh->nblocks; i++) {
			const node_location::node_id num = std::strtoll(master_mesh->meshnames[i] + 1, nullptr, 8);
			printf("%lli\n", num);
			load_locs.insert(num);
		}
		for (int i = 0; i < node_list.size(); i++) {
			node_entry_t entry;
			entry.position = positions[i];
			entry.load = bool(load_locs.find(node_list[i]) != load_locs.end());
			entry.locality_id = positions[i] * nprocs / positions.size();
			node_dir_[node_list[i]] = entry;
		}
		auto this_dir = std::move(node_dir_);
		for (int i = 0; i < nprocs; i++) {
			printf("Sending LOAD OPEN to %i\n", i);
			futs.push_back(hpx::async<load_open_action>(opts().all_localities[i], fname, this_dir));
		}
		GET(hpx::threads::run_as_os_thread(DBFreeMultimesh, master_mesh));
		for (auto& f : futs) {
			GET(f);
		}
	} else {
		std::cout << "Could not load " << fname;
		throw;
	}
	root_ptr->reconstruct_tree();
	node_registry::clear();
	futs.clear();
	for (int i = 0; i < nprocs; i++) {
		futs.push_back(hpx::async<load_close_action>(opts().all_localities[i]));
	}
	for (auto& f : futs) {
		GET(f);
	}
}

void node_server::reconstruct_tree() {
	for (integer ci = 0; ci < NCHILD; ci++) {
		is_refined = true;
		auto cloc = my_location.get_child(ci);
		auto iter = node_dir_.find(cloc.to_id());
		children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc).get();
	}
	current_time = output_time;
	rotational_time = output_rotation_count;
}

silo_var_t::silo_var_t(const std::string& name, std::size_t nx) :
		name_(name), data_(nx * nx * nx) {
	range_.first = +std::numeric_limits<real>::max();
	range_.second = -std::numeric_limits<real>::max();
}

double&
silo_var_t::operator()(int i) {
	return data_[i];
}

double silo_var_t::operator()(int i) const {
	return data_[i];
}

std::vector<silo_load_t> silo_load_t::decompress() {
	std::vector<silo_load_t> children;
	assert(nx > INX);
	for (int ci = 0; ci < NCHILD; ci++) {
		silo_load_t child;
		child.nx = nx / 2;
		child.vars.resize(vars.size());
		child.outflows.resize(vars.size());
		const integer xo = (ci & (1 << XDIM)) ? child.nx : 0;
		const integer yo = (ci & (1 << YDIM)) ? child.nx : 0;
		const integer zo = (ci & (1 << ZDIM)) ? child.nx : 0;
		for (int f = 0; f < vars.size(); f++) {
			child.vars[f].second.resize(child.nx * child.nx * child.nx);
			child.outflows[f].first = child.vars[f].first = vars[f].first;
			child.outflows[f].second = outflows[f].second / NCHILD;
			for (integer cx = 0; cx < child.nx; cx++) {
				for (integer cy = 0; cy < child.nx; cy++) {
					for (integer cz = 0; cz < child.nx; cz++) {
						const integer child_index = cx + child.nx * (cy + child.nx * cz);
						const integer parent_index = (cx + xo) + nx * ((cy + yo) + nx * (cz + zo));
						child.vars[f].second[child_index] = vars[f].second[parent_index];
					}
				}
			}
		}

		children.push_back(std::move(child));
	}
	vars.clear();
	outflows.clear();
	return std::move(children);
}

