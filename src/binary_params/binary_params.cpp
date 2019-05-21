/*
 * binary_params.cpp
 *
 *  Created on: Mar 8, 2019
 *      Author: dmarce1
 */

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <vector>
#include <cstring>
#include <silo.h>
#include <cmath>
#include <unordered_map>

#define NDIM 3

#define NOSTAR 0
#define STAR1 1
#define STAR2 2

template<class T>
void normalize3(T& vec) {
	double sum = 0.0;
	for (int d = 0; d < 3; d++) {
		sum += vec[d] * vec[d];
	}
	sum = 1.0 / sqrt(sum);
	for (int d = 0; d < 3; d++) {
		vec[d] *= sum;
	}
}

template<class T = double>
using array_type = std::vector<T>;

using map_type = std::unordered_map<std::string, array_type<>>;

map_type var_map_;
array_type<std::array<double, NDIM>> X_;
array_type<> dx_;
array_type<bool> in_plane_;
array_type<bool> in_loc_;
array_type<double> loc_x_;
double rho_mid_;
array_type<int> in_star_;

static DBfile* db_;

std::array<double, NDIM> find_com() {
	std::array<double, NDIM> com = { { 0, 0, 0 } };
	double mtot;
	const auto& rho = var_map_["rho"];
	for (int i = 0; i < rho.size(); i++) {
		const double V = std::pow(dx_[i], 3);
		const double M = V * rho[i];
		mtot += M;
		for (int d = 0; d < NDIM; d++) {
			com[d] += X_[i][d] * M;
		}
	}
	for (int d = 0; d < NDIM; d++) {
		com[d] /= mtot;
	}
	return com;
}

std::array<double, NDIM> find_vcom() {
	std::array<double, NDIM> vcom = { { 0, 0, 0 } };
	double mtot;
	const auto& rho = var_map_["rho"];
	const auto& sx = var_map_["sx"];
	const auto& sy = var_map_["sy"];
	const auto& sz = var_map_["sz"];
	for (int i = 0; i < rho.size(); i++) {
		const double V = std::pow(dx_[i], 3);
		const double M = V * rho[i];
		mtot += M;
		vcom[0] += sx[i] * V;
		vcom[1] += sy[i] * V;
		vcom[2] += sz[i] * V;
	}
	for (int d = 0; d < NDIM; d++) {
		vcom[d] /= mtot;
	}
	return vcom;
}

double find_omega(std::array<double, NDIM> com) {
	const auto& rho = var_map_["rho"];
	const auto& sx = var_map_["sx"];
	const auto& sy = var_map_["sy"];
	double I = 0.0;
	double l = 0.0;
	for (int i = 0; i < rho.size(); i++) {
		if (rho[i] > 1.) {
			const double V = std::pow(dx_[i], 3);
			std::array<double, NDIM> X;
			for (int d = 0; d < NDIM; d++) {
				X[d] = X_[i][d] - com[d];
			}
			const double r2 = X[0] * X[0] + X[1] * X[1];
			const double this_omega = (X[0] * sy[i] - X[1] * sx[i]) / r2
					/ rho[i];
			I += rho[i] * V * r2;
			l += (X[0] * sy[i] - X[1] * sx[i]) * V;
		}

	}
	printf("%e %e\n", l, I);
	return l / I;
}

void tag_loc(std::array<double, 2> loc, std::array<double, 3> a) {
	const auto& rho = var_map_["rho"];
	double p[3], n[3];
	n[0] = loc[0];
	n[1] = loc[1];
	n[2] = 0;
	for (int i = 0; i < rho.size(); i++) {
		const auto p = X_[i];
		double amp[3];
		double ampdotn;
		for (int d = 0; d < 3; d++) {
			amp[d] = a[d] - p[d];
		}
		ampdotn = 0.0;
		for (int d = 0; d < 3; d++) {
			ampdotn += amp[d] * n[d];
		}
		double dist = 0.0;
		double loc_x = 0.0;
		for (int d = 0; d < 3; d++) {
			dist += std::pow(amp[d] - ampdotn * n[d], 2);
			loc_x += n[d] * p[d];
		}
		dist = std::sqrt(dist);
		loc_x_.push_back(loc_x);
		if (dist < dx_[i]) {
			in_loc_.push_back(true);
		} else {
			in_loc_.push_back(false);
		}
	}
}

double find_eigenvector(std::array<double, 2>& e) {
	std::array<std::array<double, 2>, 2> q = { { { 0, 0 }, { 0, 0 } } };

	const auto& rho = var_map_["rho"];

	for (int i = 0; i < rho.size(); i++) {
		const double V = std::pow(dx_[i], 3);
		const double M = V * rho[i];
		const double x = X_[i][0];
		const double y = X_[i][1];
		q[0][0] += x * x * M;
		q[0][1] += x * y * M;
		q[1][0] += y * x * M;
		q[1][1] += y * y * M;
	}

	std::array<double, 2> b0, b1;
	double A, bdif;
	int iter = 0;
	b0[0] = 1.0;
	b0[1] = 1.0;
	do {
		iter++;
		b1[0] = b1[1] = 0.0;
		for (int i = 0; i < 2; i++) {
			for (int m = 0; m < 2; m++) {
				b1[i] += q[i][m] * b0[m];
			}
		}
		A = sqrt(b1[0] * b1[0] + b1[1] * b1[1]);
		bdif = 0.0;
		for (int i = 0; i < 2; i++) {
			b1[i] = b1[i] / A;
			bdif += pow(b0[i] - b1[i], 2);
		}
		for (int i = 0; i < 2; i++) {
			b0[i] = b1[i];
		}
	} while (fabs(bdif) > 1.0e-14);
	double lambda = 0.0;
	double e2 = 0.0;
	e = b0;
	for (int m = 0; m < 2; m++) {
		lambda += e[m] * (q[m][0] * e[0] + q[m][1] * e[1]);
		e2 += e[m] * e[m];
	}
	return lambda / e2;
}

double sum_all(const std::string var_name) {
	const auto var = var_map_[var_name];
	double sum = 0.0;
	for (int i = 0; i < var.size(); i++) {
		const auto V = dx_[i] * dx_[i] * dx_[i];
		sum += var[i] * V;
	}
	return sum;
}

std::pair<double, int> max_all(const std::string var_name, bool plane_only =
		false) {
	const auto var = var_map_[var_name];
	double max = -std::numeric_limits<double>::max();
	int max_i = 0;
	for (int i = 0; i < var.size(); i++) {
		if (!plane_only || in_plane_[i]) {
			if (max < var[i]) {
				max = var[i];
				max_i = i;
			}
		}
	}
	return std::pair<double, int>(max, max_i);
}

std::pair<double, int> min_all(const std::string var_name, bool plane_only =
		false) {
	const auto var = var_map_[var_name];
	double min = +std::numeric_limits<double>::max();
	int min_i = 0;
	for (int i = 0; i < var.size(); i++) {
		if (!plane_only || in_plane_[i]) {
			if (min > var[i]) {
				min = var[i];
				min_i = i;
			}
		}
	}
	return std::pair<double, int>(min, min_i);
}

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
	std::string filename = argv[1];
	db_ = DBOpenReal(filename.c_str(), DB_HDF5, DB_READ);

	if (db_ == nullptr) {
		printf("Unable to open %s\n", filename.c_str());
		return -1;
	}

	long long n_species;
	long long version;
	double omega;
	double code_to_s;
	double cgs_time;
	DBReadVar(db_, "cgs_time", (void*) &cgs_time);
	DBReadVar(db_, "version", (void*) &version);
	DBReadVar(db_, "n_species", (void*) &n_species);
	DBReadVar(db_, "code_to_s", (void*) &code_to_s);
	DBReadVar(db_, "omega", (void*) &omega);
	printf("Omega = %e\n", omega);
	printf("SILO version: %i\n", version);
	printf("N species   : %i\n", n_species);

	printf("Reading table of contents\n");
	DBmultimesh* mmesh = DBGetMultimesh(db_, "quadmesh");
	std::vector<std::string> dir_names;
	for (int i = 0; i < mmesh->nblocks; i++) {
		const std::string dir = strip_nonnumeric(mmesh->meshnames[i]);
		dir_names.push_back(dir);
	}
	DBFreeMultimesh(mmesh);

	for (int i = 0; i < dir_names.size(); i++) {
		const std::string dir = dir_names[i];
		if (dir == "Decomposition") {
			continue;
		}
//		printf("%i of %i - %s\n", i, dir_names.size(), dir.c_str());
		DBSetDir(db_, dir.c_str());

		int sz;

		const DBtoc* this_toc = DBGetToc(db_);
		for (int j = 0; j < this_toc->nvar; j++) {
			const std::string qvar = this_toc->qvar_names[j];
			DBquadvar* var = DBGetQuadvar(db_, qvar.c_str());
			sz = var->nels;
			auto& data = var_map_[qvar];
			data.resize(data.size() + sz);
			double* dest = &(data[data.size() - sz]);
			if (version == 100
					&& (qvar == "sx" || qvar == "sy" || qvar == "sz")) {
				for (int k = 0; k < sz; k++) {
					(((double**) var->vals)[0])[k] *= code_to_s;
				}
			}
			std::memcpy(dest, ((double**) var->vals)[0], sizeof(double) * sz);
			DBFreeQuadvar(var);
		}

		DBquadmesh* mesh = DBGetQuadmesh(db_, "quadmesh");
		const double* xc = static_cast<double*>(mesh->coords[0]);
		const double* yc = static_cast<double*>(mesh->coords[1]);
		const double* zc = static_cast<double*>(mesh->coords[2]);
		const double dx = xc[1] - xc[0];
		for (int l = 0; l < mesh->dims[2] - 1; l++) {
			for (int k = 0; k < mesh->dims[1] - 1; k++) {
				for (int j = 0; j < mesh->dims[0] - 1; j++) {
					std::array<double, NDIM> this_X = { xc[j] + 0.5 * dx, yc[k]
							+ 0.5 * dx, zc[l] + 0.5 * dx, };
					X_.push_back(this_X);
					dx_.push_back(dx);
					const bool in_plane = std::abs(zc[l] - 0.5 * dx) < dx;
					in_plane_.push_back(in_plane);
				}
			}
		}

		DBFreeQuadmesh(mesh);

		DBSetDir(db_, "..");

	}
	auto total_size = var_map_["rho_1"].size();
	auto& rho = var_map_["rho"];
	rho.resize(total_size, 0.0);

	for (int i = 0; i < n_species; i++) {
		std::string name = std::string("rho_") + char('1' + i);
		auto& this_rho = var_map_[name];
		for (int j = 0; j < total_size; j++) {
			rho[j] += this_rho[j];
		}
	}

	printf("Mass sum = %e\n\n\n\n", sum_all("rho"));

	auto rho_max = max_all("rho").first;
	auto rho_min = min_all("rho").first;
	rho_mid_ = std::sqrt(rho_max * rho_min);

	printf("rho_max = %e\n", rho_max);
	printf("rho_mid = %e\n", rho_mid_);
	printf("rho_min = %e\n", rho_min);

	auto c1i = max_all("rho_1").second;
	auto c2i = max_all("rho_3").second;

	printf("c1 at %e %e %e\n", X_[c1i][0], X_[c1i][1], X_[c1i][2]);
	printf("c2 at %e %e %e\n", X_[c2i][0], X_[c2i][1], X_[c2i][2]);

	printf("Mass sum = %e\n", sum_all("rho_1") / 2e33);
	printf("Mass sum = %e\n", sum_all("rho_2") / 2e33);
	printf("Mass sum = %e\n", sum_all("rho_3") / 2e33);
	printf("Mass sum = %e\n", sum_all("rho_4") / 2e33);
	printf("Mass sum = %e\n", sum_all("rho_5") / 2e33);

	auto com = find_com();
	auto vcom = find_vcom();

	printf("Center of Mass = %e %e %e\n", com[0], com[1], com[2]);

	std::array<double, 2> loc;
	double q = find_eigenvector(loc);
	printf("LOC = %e  %e\n", loc[0], loc[1]);
	printf("q = %e\n", q);

	tag_loc(loc, com);
	{
		FILE* fp = fopen("loc.txt", "wt");
		auto& rho = var_map_["rho"];
		for (int i = 0; i < rho.size(); i++) {
			if (in_loc_[i]) {
				fprintf(fp, "%e %e\n", loc_x_[i], rho[i]);
			}
		}
		fclose(fp);
	}

	omega = find_omega(com);
	double period = 2.0 * M_PI / omega / 60. / 60. / 24.;
	printf("Omega = %e\n", omega);
	printf("Period = %e days\n", period);

	double l1 = -std::numeric_limits<double>::max();
	double l2 = -std::numeric_limits<double>::max();
	double l3 = -std::numeric_limits<double>::max();
	double l1_loc, l2_loc, l3_loc;

	double c1_loc = X_[c1i][0] * loc[0] + X_[c1i][1] * loc[1];
	double c2_loc = X_[c2i][0] * loc[0] + X_[c2i][1] * loc[1];

	{
		const auto& phi = var_map_["phi"];
		for (int i = 0; i < phi.size(); i++) {
			double this_loc = X_[i][0] * loc[0] + X_[i][1] * loc[1];
			if (in_loc_[i]) {
				double phi_eff = -omega
						* (X_[i][0] * X_[i][0] + X_[i][1] * X_[i][1]) + phi[i];
				double* lptr;
				double* loc_ptr;
				if ((c1_loc < c2_loc && this_loc < c1_loc)
						|| (c1_loc > c2_loc && this_loc > c1_loc)) {
					lptr = &l3;
					loc_ptr = &l3_loc;
				} else if ((this_loc - c1_loc) * (this_loc - c2_loc) <= 0.0) {
					lptr = &l1;
					loc_ptr = &l1_loc;
				} else {
					lptr = &l2;
					loc_ptr = &l2_loc;
				}
				if (phi_eff > *lptr) {
					*lptr = phi_eff;
					*loc_ptr = this_loc;
				}
			}
		}
	}

	printf("L1 = %e @ %e\n", l1, l1_loc);
	printf("L2 = %e @ %e\n", l2, l2_loc);
	printf("L3 = %e @ %e\n", l3, l3_loc);

	const auto& sx = var_map_["sx"];
	const auto& sy = var_map_["sy"];
	const auto& sz = var_map_["sz"];
	const auto& gx = var_map_["gx"];
	const auto& gy = var_map_["gy"];
	const auto& gz = var_map_["gz"];
	const auto& tau = var_map_["tau"];

#define NREGION 3

	double M[NREGION] = { 0, 0, 0 };
	double spin[NREGION] = { 0, 0, 0 };
	double lorb = 0.0;
	system( "rm star2.txt");
	FILE* fp = fopen("roche.txt", "wt");
	double T = 0.0, W = 0.0;
	{

		const auto& phi = var_map_["phi"];
		for (int i = 0; i < phi.size(); i++) {
			double R2 = X_[i][0] * X_[i][0] + X_[i][1] * X_[i][1];
			double phi_eff = phi[i] - 0.5 * R2 * omega;
			double g1 = 0, g2 = 0;
			std::array<double, 3> d1, d2;
			for (int d = 0; d < NDIM; d++) {
				d1[d] = X_[c1i][d] - X_[i][d];
				d2[d] = X_[c2i][d] - X_[i][d];
			}
			normalize3(d1);
			normalize3(d2);
			g1 += (gx[i] + omega * omega * X_[i][0]) * d1[0];
			g1 += (gy[i] + omega * omega * X_[i][1]) * d1[1];
			g1 += (gz[i]) * d1[2];
			g2 += (gx[i] + omega * omega * X_[i][0]) * d2[0];
			g2 += (gy[i] + omega * omega * X_[i][1]) * d2[1];
			g2 += (gz[i]) * d2[2];
			int v;
			std::array<double, 3> pivot;
			if (g1 >= 0 && g2 >= 0) {
				if (g1 > g2) {
					v = STAR1;
				} else {
					v = STAR2;
				}
			} else if (g1 >= 0) {
				v = STAR1;
			} else if (g2 >= 0) {
				v = STAR2;
			} else {
				v = NOSTAR;
			}
			in_star_.push_back(v);
		}

		std::array<std::array<double, NDIM>, NREGION> this_com;
		std::array<std::array<double, NDIM>, NREGION> this_vcom;

		for (int f = 0; f < NREGION; f++) {
			for (int d = 0; d < NDIM; d++) {
				this_com[f][d] = this_vcom[f][d];
			}
		}

		for (int i = 0; i < phi.size(); i++) {
			const auto v = in_star_[i];
			const double vol = dx_[i] * dx_[i] * dx_[i];
			for (int d = 0; d < NDIM; d++) {
				this_com[v][d] += rho[i] * X_[i][d] * vol;
			}
			this_vcom[v][0] += sx[i] * vol;
			this_vcom[v][1] += sy[i] * vol;
			this_vcom[v][2] += sz[i] * vol;
			M[v] += rho[i] * vol;
		}

		for (int f = 0; f < NREGION; f++) {
			for (int d = 0; d < NDIM; d++) {
				this_com[f][d] /= M[f];
				this_vcom[f][d] /= M[f];
			}
		}


		for (int i = 0; i < phi.size(); i++) {
			const auto v = in_star_[i];

			const double vol = dx_[i] * dx_[i] * dx_[i];
			std::array<double, NDIM> R;
			for (int d = 0; d < NDIM; d++) {
				R[d] = X_[i][d] - this_com[v][d];
			}
			const auto ds = (sy[i] * R[0] - sx[i] * R[1]) * vol;
			spin[v] += ds;
			if (v == STAR1 || v == STAR2) {
				for (int d = 0; d < NDIM; d++) {
					R[d] = X_[i][d] - com[d];
				}
				const auto ds = (sy[i] * R[0] - sx[i] * R[1]) * vol;
				lorb += ds;
			}

			if( v == STAR2 && rho[i] > 1.0e-6 ) {
				double R2 = 0.0;
				for (int d = 0; d < 2; d++) {
					R[d] = X_[i][d] - this_com[v][d];
					R2 += R[d] * R[d];
				}
				const auto this_omega = ((sy[i] * R[0] - sx[i] * R[1]) / rho[i] - (this_vcom[v][1] * R[0] - this_vcom[v][0] * R[1])) / R2;
				T += std::pow(this_omega, 2) * R2 * rho[i] * vol;
				W += phi[i] * rho[i] * vol;
			}

		}


	}

	lorb -= spin[STAR2] + spin[STAR1];
	fclose(fp);
#define sqr(a) ((a)*(a))

	for (int f = 0; f < NREGION; f++) {
		M[f] /= 1.99e+33;
	}
	double sep =
			std::sqrt(
					sqr(
							X_[c1i][0]
									- X_[c2i][0]) + sqr(X_[c1i][1] - X_[c2i][1]) + sqr(X_[c1i][2] - X_[c2i][2]));
	fp = fopen("binary.txt", "at");
	fprintf(fp, "%e %e %e %e %e %e %e %e %e\n", cgs_time, M[0], M[1], M[2], sep,
			spin[1], spin[2], lorb, T / -W);
	fclose(fp);
	/* Close SILO */

	DBClose(db_);

}

/*
 using space_vector = std::array<double,3>;
 using std::sqrt;

 template<class T>
 inline T sqr(T s) {
 return s * s;
 }

 std::string strip_nonnumeric(std::string&& s) {
 s.erase(std::remove_if(s.begin(), s.end(), [](char c) {return c < '0' || c > '9';}), s.end());
 return std::move(s);
 }

 double find_eigenvector(const std::array<std::array<double, 3>, 3>& q, std::array<double, 3>& e) {
 std::array<double, 3> b0, b1;
 double A, bdif;
 int iter = 0;
 b0[2] = 0.0;
 b0[0] = 1.0;
 b0[1] = 1.0;
 do {
 iter++;
 b1[0] = b1[1] = b1[2] = 0.0;
 for (int i = 0; i < 3; i++) {
 for (int m = 0; m < 3; m++) {
 b1[i] += q[i][m] * b0[m];
 }
 }
 A = sqrt(sqr(b1[0]) + sqr(b1[1]) + sqr(b1[2]));
 bdif = 0.0;
 for (int i = 0; i < 3; i++) {
 b1[i] = b1[i] / A;
 bdif += pow(b0[i] - b1[i], 2);
 }
 for (int i = 0; i < 3; i++) {
 b0[i] = b1[i];
 }
 } while (fabs(bdif) > 1.0e-14);
 double lambda = 0.0;
 double e2 = 0.0;
 e = b0;
 for (int m = 0; m < 3; m++) {
 lambda += e[m] * (q[m][0]*e[0]+q[m][1]*e[1]+q[m][2]*e[2]);
 e2 += e[m] * e[m];
 }
 return lambda / e2;
 }


 std::array<double, 3> center_of_mass(const std::vector<cell_t>& cells) {
 std::array<double, 3> com;
 double mass = 0.0;
 for (int d = 0; d < 3; d++) {
 com[d] = 0.0;
 }
 for (const auto& c : cells) {
 const auto this_vol = c.dx * c.dx * c.dx;
 double rho = 0.0;
 for (int i = 0; i < cell_t::n_species; i++) {
 rho += c.rho[i];
 }
 const auto this_mass = this_vol * rho;
 mass += this_mass;
 for (int d = 0; d < 3; d++) {
 com[d] += this_mass * c.x[d];
 }
 }
 for (int d = 0; d < 3; d++) {
 com[d] /= mass;
 }
 return com;
 }

 double total_mass(std::vector<cell_t>& cells) {
 double mass = 0.0;
 for (auto& c : cells) {
 const auto this_vol = c.dx * c.dx * c.dx;
 double rho = 0.0;
 for (int i = 0; i < cell_t::n_species; i++) {
 rho += c.rho[i];
 }
 c.rho_tot = rho;
 const auto this_mass = this_vol * rho;
 mass += this_mass;
 }
 return mass;
 }

 std::array<std::array<double, 3>, 3> quadrupole_moment(const std::vector<cell_t>& cells, const std::array<double, 3>& com) {
 std::array<std::array<double, 3>, 3> q = { { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } };
 for (const auto& c : cells) {
 for (int n = 0; n < 3; n++) {
 const auto this_vol = c.dx * c.dx * c.dx;
 double rho = 0.0;
 for (int i = 0; i < cell_t::n_species; i++) {
 rho += c.rho[i];
 }
 const auto this_mass = this_vol * rho;
 double x = c.x[0] - com[0];
 double y = c.x[1] - com[1];
 double z = c.x[2] - com[2];
 double r2 = x * x + y * y;
 q[0][0] += 3 * this_mass * x * x;
 q[0][1] += 3 * this_mass * x * y;
 //			q[0][2] += 3 * this_mass * x * z;
 q[1][0] += 3 * this_mass * y * x;
 q[1][1] += 3 * this_mass * y * y;
 //			q[1][2] += 3 * this_mass * y * z;
 //			q[2][0] += 3 * this_mass * z * x;
 //			q[2][1] += 3 * this_mass * z * y;
 //			q[2][2] += 3 * this_mass * z * z;
 q[0][0] -= r2 * this_mass;
 q[1][1] -= r2 * this_mass;
 //			q[2][2] -= r2 * this_mass;
 }
 }
 return q;
 }

 int main(int argc, char* argv[]) {

 std::set<std::string> mesh_names;
 std::set<std::string> var_names;
 std::vector<cell_t> cells;
 var_names.insert("phi");
 var_names.insert("sx");
 var_names.insert("sy");
 var_names.insert("sz");

 if (argc != 2) {
 printf("Usage: binary_params <silo_file>\n");
 return -1;
 }

 std::string filename = argv[1];
 auto handle = DBOpenReal(filename.c_str(), DB_HDF5, DB_READ);

 if (handle == nullptr) {
 printf("Unable to open %s\n", filename.c_str());
 return -1;
 }

 DBReadVar(handle, "n_species", &cell_t::n_species);
 for (int i = 0; i < cell_t::n_species; i++) {
 var_names.insert("rho_" + std::to_string(i + 1));
 }
 auto mesh = DBGetMultimesh(handle, "quadmesh");
 for (int i = 0; i < mesh->nblocks; i++) {
 auto name = strip_nonnumeric(std::string(mesh->meshnames[i]));
 mesh_names.insert(name);
 }

 bool first_call = true;
 for (const auto& vn : var_names) {
 printf("Reading %s\n", vn.c_str());
 int p = 0;
 for (auto const& mn : mesh_names) {
 std::string mesh_loc = "/" + mn + "/quadmesh";
 std::string var_loc = "/" + mn + "/" + vn;

 auto quadmesh = DBGetQuadmesh(handle, mesh_loc.c_str());

 double* X = (double*) quadmesh->coords[0];
 const auto dx = X[1] - X[0];
 const auto dv = dx * dx * dx;

 //		printf("Reading %s\n", var_loc.c_str());
 auto var = DBGetQuadvar(handle, var_loc.c_str());
 if (var == nullptr) {
 printf("Unable to read %s\n", var_loc.c_str());
 return -1;
 }
 int i = 0;
 for (int l = 0; l < var->dims[2]; l++) {
 for (int k = 0; k < var->dims[1]; k++) {
 for (int j = 0; j < var->dims[0]; j++) {
 if (first_call) {
 cell_t c;
 c.dx = dx;
 c.x[0] = ((double*) quadmesh->coords[0])[j] + 0.5 * dx;
 c.x[1] = ((double*) quadmesh->coords[1])[k] + 0.5 * dx;
 c.x[2] = ((double*) quadmesh->coords[2])[l] + 0.5 * dx;
 for (int d = 0; d < 3; d++) {
 //			printf("%e ", c.x[d]);
 }
 cells.push_back(c);
 }
 auto& c = cells[p];
 const auto& val = ((double*) var->vals[0])[i];
 //		printf("%e\n", val);
 if (vn == std::string("phi")) {
 c.phi = val;
 } else {
 if (vn.size() == 2 && vn[0] == 's') {
 int index = vn[1] - 'x';
 c.s[index] = val;
 } else if (std::strncmp("rho", vn.c_str(), 3) == 0) {
 int index = vn[4] - '1';
 c.rho[index] = val;
 } else {
 printf("Error on line %i\n", __LINE__);
 return -1;
 }
 }
 p++;
 i++;
 }
 }
 }

 DBFreeQuadvar(var);
 DBFreeQuadmesh(quadmesh);
 }
 first_call = false;
 }
 DBClose(handle);

 auto M = total_mass(cells);
 printf("Total Mass: %e\n", M);

 auto com = center_of_mass(cells);
 printf("Center of Mass: %e %e %e\n", com[0], com[1], com[2]);

 auto q = quadrupole_moment(cells, com);
 printf("Quadrupole Moment: %12e %12e\n", q[0][0], q[0][1]);
 printf("                   %12e %12e\n", q[1][0], q[1][1]);

 double lambda;
 std::array<double, 3> loc;
 lambda = find_eigenvector(q, loc);

 printf("Line of Centers:   %12e %12e %12e\n", loc[0], loc[1], loc[2]);

 double rho_max = 0.0;
 space_vector c1, c2;
 for( const auto& c : cells) {
 double dx2_max = 0.0;
 for( int d = 0; d < 3; d++) {
 dx2_max = std::max(dx2_max,sqr(c.x[d] - loc[d]));
 }
 auto rho = c.rho_tot;
 if( dx2_max <= c.dx* c.dx ) {
 if( rho > rho_max) {
 rho_max = rho;
 c1 = c.x;
 }
 }
 }


 double d1 = sqrt(c1[0]*c1[0]+c1[1]*c1[1]+c1[2]*c1[2]);
 double d2 = lambda / d1 / M;

 double a = d1 + d2;
 printf( "%e %e %e\n", d1, d2, a);

 for( int d= 0; d < 3; d++) {
 c2[d] = c1[d] + loc[d]*(d1+d2);
 }
 printf( "First  star at %e %e %e with rho_max = %e\n", c1[0], c1[1], c1[2], rho_max);
 printf( "Second star at %e %e %e\n", c2[0], c2[1], c2[2]);


 return 0;
 }
 */
