/*
 * silo_planes.cpp
 *
 *  Created on: May 16, 2019
 *      Author: dmarce1
 */

#include <silo.h>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>

#define SILO_DRIVER DB_HDF5
#define ZERO (1e-10)

template<class T = double>
using array_type = std::vector<T>;

using map_type = std::unordered_map<std::string, array_type<>>;

map_type var_map_;
array_type<bool> in_plane_;
array_type<bool> in_loc_;
array_type<double> loc_x_;
double rho_mid_;
array_type<int> in_star_;

DBfile* db_in_;
DBfile* db_out_;

std::string strip_nonnumeric(std::string&& s) {
	s.erase(
			std::remove_if(s.begin(), s.end(),
					[](char c) {return c < '0' || c > '9';}), s.end());
	return std::move(s);
}

int main(int argc, char* argv[]) {

	if (argc != 2) {
		printf("Usage: binary_params <silo_file>\n");
		return -1;
	}

	/* Open SILO */

	printf("Opening SILO\n");
	std::string in_filename = argv[1];
	std::string out_filename = std::string("plane.") + in_filename;

	db_in_ = DBOpenReal(in_filename.c_str(), DB_HDF5, DB_READ);
	db_out_ = DBCreateReal(out_filename.c_str(), DB_CLOBBER, DB_LOCAL,
			"Octo-tiger", SILO_DRIVER);

	if (db_in_ == nullptr) {
		printf("Unable to open %s\n", in_filename.c_str());
		return -1;
	}

	long long n_species;
	long long version;
	double omega;
	double code_to_s;
	double cgs_time;
	DBReadVar(db_in_, "cgs_time", (void*) &cgs_time);
	DBReadVar(db_in_, "version", (void*) &version);
	DBReadVar(db_in_, "n_species", (void*) &n_species);
	DBReadVar(db_in_, "code_to_s", (void*) &code_to_s);
	DBReadVar(db_in_, "omega", (void*) &omega);
	printf("Omega = %e\n", omega);
	printf("SILO version: %i\n", version);
	printf("N species   : %i\n", n_species);

	printf("Reading table of contents\n");

	DBmultimesh* mmesh = DBGetMultimesh(db_in_, "quadmesh");
	std::vector<std::string> dir_names;
	for (int i = 0; i < mmesh->nblocks; i++) {
		const std::string dir = strip_nonnumeric(mmesh->meshnames[i]);
		dir_names.push_back(dir);
	}
	DBFreeMultimesh(mmesh);

	int n_total_domains = 0;
	std::vector<char*> mesh_names;
	std::vector<std::vector<char*>> var_names;

	std::vector<std::string> top_var_names;

	bool first_pass = true;



	const DBtoc* toc = DBGetToc(db_in_);

	for( int i = 0; i < toc->nvar; i++) {
		int one = 1;
		const auto vlen = DBGetVarByteLength (db_in_, toc->var_names[i]);
		const auto type = DBGetVarType(db_in_, toc->var_names[i]);
		void* data = malloc(vlen);
		DBReadVar(db_in_, toc->var_names[i], data);
		DBWrite(db_out_, toc->var_names[i], data, &one, 1, type);
		free(data);
	}

	for (int i = 0; i < dir_names.size(); i++) {
		const std::string dir = dir_names[i];
		if (dir == "Decomposition") {
			continue;
		}
		DBSetDir(db_in_, dir.c_str());

		int sz;

		const DBtoc* this_toc = DBGetToc(db_in_);


		DBquadmesh* mesh = DBGetQuadmesh(db_in_, "quadmesh");
		const double* xc = static_cast<double*>(mesh->coords[0]);
		const double* yc = static_cast<double*>(mesh->coords[1]);
		const double* zc = static_cast<double*>(mesh->coords[2]);
		const double dx = xc[1] - xc[0];
		int l_plane = -1;

		for (int l = 0; l < mesh->dims[2] - 1; l++) {
			if (zc[l] > -0.5*dx && zc[l] <= 0.5*dx) {
				l_plane = 1;
				break;
			}
		}

		if (l_plane != -1) {
			DBMkDir(db_out_, dir.c_str());
			DBSetDir(db_out_, dir.c_str());
			int one = 1;
			auto optlist_var = DBMakeOptlist(1);
			DBAddOption(optlist_var, DBOPT_HIDE_FROM_GUI, &one);

			const auto name = "quadmesh_2d";
			const auto& coord_names = mesh->labels;
			const auto& coords = mesh->coords;
			const auto& dims = mesh->dims;
			const int ndims = 2;
			const auto datatype = mesh->datatype;
			const auto coordtype = mesh->coordtype;
			DBPutQuadmesh(db_out_, "quadmesh_2d", mesh->labels, mesh->coords,
					mesh->dims, 2, mesh->datatype, mesh->coordtype,
					optlist_var);
			std::string mesh_name = dir + "/quadmesh_2d";

			char* str = new char[mesh_name.size() + 1];
			strcpy(str, mesh_name.c_str());
			mesh_names.push_back(str);
			n_total_domains++;

			for (int j = 0; j < this_toc->nqvar; j++) {
				const std::string qvar = this_toc->qvar_names[j];
				DBquadvar* var = DBGetQuadvar(db_in_, qvar.c_str());
				if (first_pass) {
					top_var_names.push_back(var->name);
					var_names.push_back(std::vector<char*>());
				}
				DBPutQuadvar1(db_out_, var->name, "quadmesh_2d", var->vals[0],
						var->dims, 2, (const void*) NULL, 0, DB_DOUBLE,
						DB_ZONECENT, optlist_var);

				std::string var_name = dir + "/" + var->name;
				char* str = new char[var_name.size() + 1];
				strcpy(str, var_name.c_str());

				var_names[j].push_back(str);

				DBFreeQuadvar(var);

			}
			first_pass = false;

			DBFreeOptlist(optlist_var);
			DBSetDir(db_out_, "..");

		}

		DBFreeQuadmesh(mesh);

		DBSetDir(db_in_, "..");

	}

	int mesh_type = DB_QUADMESH;
	auto optlist_var = DBMakeOptlist(1);
	DBAddOption(optlist_var, DBOPT_MB_BLOCK_TYPE, &mesh_type);
	DBPutMultimesh(db_out_, "quadmesh_2d", n_total_domains, mesh_names.data(),
	NULL, optlist_var);
	DBFreeOptlist(optlist_var);
	for (int f = 0; f < top_var_names.size(); f++) {
		DBPutMultivar(db_out_, top_var_names[f].c_str(), n_total_domains,
				var_names[f].data(),
				std::vector<int>(n_total_domains, DB_QUADVAR).data(), NULL);
	}


	auto tmp = DBGetDefvars(db_in_, "expressions");
	DBPutDefvars(db_out_, "expressions", tmp->ndefs, tmp->names, tmp->types, tmp->defns, NULL);
	DBFreeDefvars(tmp);


	DBClose(db_in_);
	DBClose(db_out_);
	return 0;
}
