#include "octotiger/unitiger/hydro_impl/reconstruct_kernel_interface.hpp"

#pragma GCC push_options
#pragma GCC optimize("unroll-loops")

#include <Vc/Vc>
#include <Vc/common/mask.h>
#include <Vc/vector.h>

using vc_type = Vc::Vector<double, Vc::VectorAbi::Avx>;
using mask_type = vc_type::mask_type;
using index_type = Vc::Vector<int, Vc::VectorAbi::Avx>;

template <>
inline vc_type copysign_wrapper<vc_type>(const vc_type& tmp1, const vc_type& tmp2) {
    return Vc::copysign(tmp1, tmp2);
}
template <>
inline vc_type abs_wrapper<vc_type>(const vc_type& tmp1) {
    return Vc::abs(tmp1);
}
template <>
inline vc_type minmod_wrapper<vc_type>(const vc_type& a, const vc_type& b) {
    return (copysign_wrapper<vc_type>(0.5, a) + copysign_wrapper<vc_type>(0.5, b)) *
        min_wrapper<vc_type>(abs_wrapper<vc_type>(a), abs_wrapper<vc_type>(b));
}
template <>
inline vc_type minmod_theta_wrapper<vc_type>(const vc_type& a, const vc_type& b, const vc_type& c) {
    return minmod_wrapper<vc_type>(c * minmod_wrapper<vc_type>(a, b), 0.5 * (a + b));
}
template <>
inline vc_type load_value<vc_type>(const double* __restrict__ data, const size_t index) {
    return vc_type(data + index);
}

inline void make_monotone_wrapper(double* __restrict__ qld, const double* __restrict__ q0d,
    double* __restrict__ qrd, const mask_type& mask) {
    vc_type ql(qld);
    vc_type qr(qrd);
    const vc_type q0(q0d);
    const vc_type tmp1 = qr - ql;
    const vc_type tmp2 = qr + ql;

    const mask_type mask1_tmp1 = mask_type(qr < q0);
    const mask_type mask1_tmp2 = mask_type(q0 < ql);
    const mask_type mask1 = (mask1_tmp1 ^ mask1_tmp2);
    const vc_type tmp3 = tmp1 * tmp1 / 6.0;
    const vc_type tmp4 = tmp1 * (q0 - 0.5 * tmp2);
    const mask_type mask2 = mask_type(tmp4 > tmp3);
    const mask_type mask3 = !mask_type(tmp4 > tmp3) && mask_type(-tmp3 > tmp4);
    Vc::where(mask2, ql) = (3.0 * q0 - 2.0 * qr);
    Vc::where(mask3, qr) = (3.0 * q0 - 2.0 * ql);
    Vc::where(mask1, qr) = q0;
    Vc::where(mask1, ql) = q0;
    ql.store(qld);
    qr.store(qrd);
}

constexpr int number_faces = 15;
constexpr int number_dirs = 27;
constexpr int q_inx = INX + 2;
constexpr int q_inx3 = q_inx * q_inx * q_inx;
constexpr int q_face_offset = number_dirs * q_inx3;
constexpr int q_dir_offset = q_inx3;

int to_q_index(int j, int k, int l) {
    return j * q_inx * q_inx + k * q_inx + l;
}
// template <>
// inline void store_value<vc_type>(const double* __restrict__ data, const size_t index, const
// vc_type &value) {
//    return value.store(data + index);
//}

void reconstruct_minmod(double* __restrict__ q, const std::vector<safe_real>& u, int d, int f) {
    static const cell_geometry<NDIM, INX> geo;
    static constexpr auto dir = geo.direction();
    const auto di = dir[d];
    const int start_index = f * q_face_offset + d * q_dir_offset;
    for (int j = 0; j < geo.H_NX_XM4; j++) {
        for (int k = 0; k < geo.H_NX_YM4; k++) {
            for (int l = 0; l < geo.H_NX_ZM4; l++) {
                const int i = geo.to_index(j + 2, k + 2, l + 2);
                const int q_i = to_q_index(j, k, l) + start_index;
                q[q_i] = u[i] + 0.5 * minmod(u[i + di] - u[i], u[i] - u[i - di]);
            }
        }
    }
}

void reconstruct_ppm_experimental(double* __restrict__ q,
    const std::vector<safe_real>& u, bool smooth, bool disc_detect,
    const std::vector<std::vector<double>>& disc, int d, int f) {
    static const cell_geometry<NDIM, INX> geo;
    static constexpr auto dir = geo.direction();
    // const vc_type zindices = vc_type::IndexesFromZero() + 1;
    static thread_local auto D1 = std::vector<safe_real>(geo.H_N3, 0.0);
    const auto di = dir[d];
    const auto flipped_di = geo.flip(d);

    for (int j = 0; j < geo.H_NX_XM2; j++) {
        for (int k = 0; k < geo.H_NX_YM2; k++) {
            for (int l = 0; l < geo.H_NX_ZM2; l += vc_type::size()) {
                const int i = geo.to_index(j + 1, k + 1, l + 1);
                const vc_type u_plus_di(u.data() + i + di);
                const vc_type u_zero(u.data() + i);
                const vc_type diff_u = u_plus_di - u_zero;
                const vc_type u_minus_di(u.data() + i - di);
                const vc_type d1 = minmod_theta_wrapper<vc_type>(diff_u, u_zero - u_minus_di, 2.0);
                d1.store(D1.data() + i);
            }
        }
    }
    // TODO How to handle the flipped stuff? I need a di wrapper yet again
    // TODO Basically I need a q_di and a q_flipped d
    // TODO Wait, do I even need this! q is stable after all
    const int start_index = f * q_face_offset + d * q_dir_offset;
    const int start_index_flipped = f * q_face_offset + flipped_di * q_dir_offset;
    for (int j = 0; j < geo.H_NX_XM4; j++) {
        for (int k = 0; k < geo.H_NX_YM4; k++) {
            for (int l = 0; l < geo.H_NX_ZM4; l += vc_type::size()) {
                const int i = geo.to_index(j + 2, k + 2, l + 2);
                const int q_i = to_q_index(j, k, l);
                const vc_type u_plus_di(u.data() + i + di);
                const vc_type u_zero(u.data() + i);
                const vc_type diff_u = u_plus_di - u_zero;

                const vc_type u_minus_di(u.data() + i - di);
                const vc_type d1(D1.data() + i);
                const vc_type d1_plus(D1.data() + i + di);
                const vc_type d1_minus(D1.data() + i - di);

                const vc_type results = 0.5 * (u_zero + u_plus_di) + (1.0 / 6.0) * (d1 - d1_plus);
                const vc_type results_flipped =
                    0.5 * (u_minus_di + u_zero) + (1.0 / 6.0) * (d1_minus - d1);

                results.store(q + start_index + q_i);
                results_flipped.store(q + start_index_flipped + q_i);
            }
        }
    }

    /*if (experiment == 1) {
        for (int j = 0; j < geo.H_NX_XM2; j++) {
            for (int k = 0; k < geo.H_NX_YM2; k++) {
#pragma ivdep
                for (int l = 0; l < geo.H_NX_ZM2; l++) {
                    const int i = geo.to_index(j + 1, k + 1, l + 1);
                    for (int gi = 0; gi < geo.group_count(); gi++) {
                        safe_real sum = 0.0;
                        for (int n = 0; n < geo.group_size(gi); n++) {
                            const auto pair = geo.group_pair(gi, n);
                            sum += q[pair.second][i + pair.first];
                        }
                        sum /= safe_real(geo.group_size(gi));
                        for (int n = 0; n < geo.group_size(gi); n++) {
                            const auto pair = geo.group_pair(gi, n);
                            q[pair.second][i + pair.first] = sum;
                        }
                    }
                }
            }
        }
        for (int d = 0; d < geo.NDIR; d++) {
            const auto di = dir[d];
            for (int j = 0; j < geo.H_NX_XM2; j++) {
                for (int k = 0; k < geo.H_NX_YM2; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM2; l++) {
                        const int i = geo.to_index(j + 1, k + 1, l + 1);
                        const auto mx = std::max(u[i + di], u[i]);
                        const auto mn = std::min(u[i + di], u[i]);
                        q[d][i] = std::min(mx, q[d][i]);
                        q[d][i] = std::max(mn, q[d][i]);
                    }
                }
            }
        }
    } */
    if (disc_detect) {
        constexpr auto eps = 0.01;
        constexpr auto eps2 = 0.001;
        constexpr auto eta1 = 20.0;
        constexpr auto eta2 = 0.05;
        const auto di = dir[d];
        for (int j = 0; j < geo.H_NX_XM4; j++) {
            for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                for (int l = 0; l < geo.H_NX_ZM4; l++) {
                    const int i = geo.to_index(j + 2, k + 2, l + 2);
                    const auto& up = u[i + di];
                    const auto& u0 = u[i];
                    const auto& um = u[i - di];
                    const auto dif = up - um;
                    if (std::abs(dif) > disc[d][i] * std::min(std::abs(up), std::abs(um))) {
                        if (std::min(std::abs(up), std::abs(um)) /
                                std::max(std::abs(up), std::abs(um)) >
                            eps2) {
                            const auto d2p = (1.0 / 6.0) * (u[i + 2 * di] + u0 - 2.0 * u[i + di]);
                            const auto d2m = (1.0 / 6.0) * (u0 + u[i - 2 * di] - 2.0 * u[i - di]);
                            if (d2p * d2m < 0.0) {
                                double eta = 0.0;
                                if (std::abs(dif) > eps * std::min(std::abs(up), std::abs(um))) {
                                    eta = -(d2p - d2m) / dif;
                                }
                                eta = std::max(0.0, std::min(eta1 * (eta - eta2), 1.0));
                                if (eta > 0.0) {
                                    const int q_i = to_q_index(j, k, l);
                                    auto ul =
                                        um + 0.5 * minmod_theta(u[i] - um, um - u[i - 2 * di], 2.0);
                                    auto ur =
                                        up - 0.5 * minmod_theta(u[i + 2 * di] - up, up - u[i], 2.0);
                                    // auto& qp = q[d][i];
                                    // auto& qm = q[geo.flip(d)][i];
                                    auto& qp = q[start_index + q_i];
                                    auto& qm = q[start_index_flipped + q_i];
	                									qp += eta * (ur - qp);
									                 	qm += eta * (ul - qm);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (!smooth) {
        const vc_type zindices = vc_type::IndexesFromZero();
        for (int j = 0; j < geo.H_NX_XM4; j++) {
            for (int k = 0; k < geo.H_NX_YM4; k++) {
                for (int l = 0; l < geo.H_NX_ZM4; l += vc_type::size()) {
                    const int i = geo.to_index(j + 2, k + 2, l + 2);
                    const int q_i = to_q_index(j, k, l);
                    const mask_type mask = (zindices + l < geo.H_NX_ZM4);
                    make_monotone_wrapper(
                        q + start_index + q_i, u.data() + i, q + start_index_flipped + q_i, mask);
                    // auto& qp = q[geo.flip(d)][i];
                    // auto& qm = q[d][i];
                    // make_monotone(qm, u[i], qp);
                }
            }
        }
    }
}

void reconstruct_experimental(const hydro::state_type& U_,
    const hydro::x_type& X, safe_real omega, const size_t nf_, const int angmom_index_,
    const std::vector<bool>& smooth_field_, const std::vector<bool>& disc_detect_, double *__restrict__ combined_q) {
    static const cell_geometry<NDIM, INX> geo;
    static thread_local std::vector<std::vector<safe_real>> AM(
        geo.NANGMOM, std::vector<safe_real>(geo.H_N3));
    //static thread_local std::vector<std::vector<std::vector<safe_real>>> Q(
    //    nf_, std::vector<std::vector<safe_real>>(geo.NDIR, std::vector<safe_real>(geo.H_N3)));

    //std::vector<double, recycler::recycle_allocator_cuda_host<double>> combined_q(
    //    15 * 27 * 10 * 10 * 10 + 32);

    static constexpr auto xloc = geo.xloc();
    static constexpr auto levi_civita = geo.levi_civita();
    static constexpr auto vw = geo.volume_weight();
    static constexpr auto dir = geo.direction();
    const int n_species_ = physics<NDIM>::get_n_species();

    const auto dx = X[0][geo.H_DNX] - X[0][0];
    const auto& U = physics<NDIM>::pre_recon<INX>(U_, X, omega, angmom_index_ != -1);
    const auto& cdiscs = physics<NDIM>::find_contact_discs<INX>(U_);
    assert(angmom_index_ > -1);
    assert(NDIM > -1);
    const int sx_i = angmom_index_;
    const int zx_i = sx_i + NDIM;

    for (int n = 0; n < geo.NANGMOM; n++) {
        for (int j = 0; j < geo.H_NX_XM4; j++) {
            for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                for (int l = 0; l < geo.H_NX_ZM4; l++) {
                    const int i = geo.to_index(j + 2, k + 2, l + 2);
                    AM[n][i] = U[zx_i + n][i] * U[0][i];
                }
            }
        }
    }

    for (int d = 0; d < geo.NDIR; d++) {
        if (d < geo.NDIR / 2) {
            for (int f = 0; f < angmom_index_; f++) {
                reconstruct_ppm_experimental(
                    combined_q, U[f], smooth_field_[f], disc_detect_[f], cdiscs, d, f);
            }
            for (int f = sx_i; f < sx_i + NDIM; f++) {
                reconstruct_ppm_experimental(combined_q, U[f], true, false, cdiscs, d, f);
            }
        }
        for (int f = zx_i; f < zx_i + geo.NANGMOM; f++) {
            reconstruct_minmod(combined_q, U[f], d, f);
        }
        if (d < geo.NDIR / 2) {
            for (int f = angmom_index_ + geo.NANGMOM + NDIM; f < nf_; f++) {
                reconstruct_ppm_experimental(
                    combined_q, U[f], smooth_field_[f], disc_detect_[f], cdiscs, d, f);
            }
        }

        if (d != geo.NDIR / 2) {
            const int start_index_sx = sx_i * q_face_offset + d * q_dir_offset;
            const int start_index_rho = d * q_dir_offset;
            for (int n = 0; n < geo.NANGMOM; n++) {
                const vc_type zindices = vc_type::IndexesFromZero();
                for (int m = 0; m < NDIM; m++) {
                    for (int q = 0; q < NDIM; q++) {
                        const auto lc = levi_civita[n][m][q];
                        if (lc != 0) {
                            for (int j = 0; j < geo.H_NX_XM4; j++) {
                                for (int k = 0; k < geo.H_NX_YM4; k++) {
                                    for (int l = 0; l < geo.H_NX_ZM4; l += vc_type::size()) {
                                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                                        const int q_i = to_q_index(j, k, l);
                                        const vc_type results = vc_type(AM[n].data() + i) -
                                            vw[d] * lc * 0.5 * xloc[d][m] *
                                                //vc_type(Q[sx_i + q][d].data() + i) *
                                                //vc_type(Q[0][d].data() + i) * dx;
                                                vc_type(combined_q + start_index_sx + q_i) *
                                                vc_type(combined_q + start_index_rho + q_i) * dx;
                                        results.store(AM[n].data() + i);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int d = 0; d < geo.NDIR / 2; d++) {
        const auto di = dir[d];

        for (int q = 0; q < NDIM; q++) {
            const auto f = sx_i + q;
						const int start_index_f = f * q_face_offset + d * q_dir_offset;
						const int start_index_flipped = f * q_face_offset + geo.flip(d) * q_dir_offset;
						const int start_index_zero = 0 * q_face_offset + d * q_dir_offset;
						const int start_index_zero_flipped = 0 * q_face_offset + geo.flip(d) * q_dir_offset;
            for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM4; l++) {
                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                        const int q_i = to_q_index(j, k, l);
                        //const auto& rho_r = Q[0][d][i];
                        //const auto& rho_l = Q[0][geo.flip(d)][i];
                        const auto& rho_r = combined_q[start_index_zero + q_i];
                        const auto& rho_l = combined_q[start_index_zero_flipped + q_i];
                        auto& qr = combined_q[start_index_f + q_i];
                        auto& ql = combined_q[start_index_flipped + q_i];
                        //auto& qr = Q[f][d][i];
                        //auto& ql = Q[f][geo.flip(d)][i];
                        const auto& ur = U[f][i + di];
                        const auto& u0 = U[f][i];
                        const auto& ul = U[f][i - di];
                        const auto b0 = qr - ql;
                        auto b = b0;
                        for (int n = 0; n < geo.NANGMOM; n++) {
                            for (int m = 0; m < NDIM; m++) {
                                const auto lc = levi_civita[n][m][q];
                                b += 12.0 * AM[n][i] * lc * xloc[d][m] / (dx * (rho_l + rho_r));
                            }
                        }
                        double blim;
                        if ((ur - u0) * (u0 - ul) <= 0.0) {
                            blim = 0.0;
                        } else {
                            blim = b0;
                        }
                        b = minmod(blim, b);
                        qr += 0.5 * (b - b0);
                        ql -= 0.5 * (b - b0);
                        if (ur > u0 && u0 > ul) {
                            if (qr > ur) {
                                ql -= (qr - ur);
                                qr = ur;
                            } else if (ql < ul) {
                                qr -= (ql - ul);
                                ql = ul;
                            }
                        } else if (ur < u0 && u0 < ul) {
                            if (qr < ur) {
                                ql -= (qr - ur);
                                qr = ur;
                            } else if (ql > ul) {
                                qr -= (ql - ul);
                                ql = ul;
                            }
                        }
                        make_monotone(qr, u0, ql);
                    }
                }
            }
        }
        if (d != geo.NDIR / 2) {
						const int start_index_sx = sx_i * q_face_offset + d * q_dir_offset;
						const int start_index_sy = sy_i * q_face_offset + d * q_dir_offset;
            for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM4; l++) {
                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                        const int q_i = to_q_index(j, k, l);
                        combined_q[start_index_sx + q_i] -= omega * (X[1][i] + 0.5 * xloc[d][1] * dx);
                        combined_q[start_index_sy + q_i] += omega * (X[0][i] + 0.5 * xloc[d][0] * dx);
                        //Q[sx_i][d][i] -= omega * (X[1][i] + 0.5 * xloc[d][1] * dx);
                        //Q[sy_i][d][i] += omega * (X[0][i] + 0.5 * xloc[d][0] * dx);
                    }
                }
            }
        }
    }
    //    for (int d = 0; d < geo.NDIR; d++) {
    // }

    // post_recon_experimental(Q, X, omega, angmom_index_ != -1);

    for (int d = 0; d < geo.NDIR; d++) {
				const int start_index_rho = rho_i * q_face_offset + d * q_dir_offset;
				const int start_index_egas = egas_i * q_face_offset + d * q_dir_offset;
				const int start_index_pot = pot_i * q_face_offset + d * q_dir_offset;
        if (d != geo.NDIR / 2) {
            for (int n = 0; n < geo.NANGMOM; n++) {
								const int start_index_sx = sx_i * q_face_offset + d * q_dir_offset;
								const int start_index_lx_n = (lx_i + n) * q_face_offset + d * q_dir_offset;
                for (int q = 0; q < NDIM; q++) {
                    for (int m = 0; m < NDIM; m++) {
                        const auto lc = levi_civita[n][m][q];
                        if (lc != 0) {
                            const vc_type xloc_tmp = vc_type(0.5 * xloc[d][m] * dx);
                            for (int j = 0; j < geo.H_NX_XM4; j++) {
                                for (int k = 0; k < geo.H_NX_YM4; k++) {
                                    for (int l = 0; l < geo.H_NX_ZM4; l += vc_type::size()) {
                                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                                        const int q_i = to_q_index(j, k, l);
                                        //const vc_type rho(Q[rho_i][d].data() + i);
                                        //const vc_type q_lx_val(Q[lx_i + n][d].data() + i);
                                        const vc_type rho(combined_q + start_index_rho + q_i);
                                        const vc_type q_lx_val(combined_q + start_index_lx_n + q_i);
                                       /* const auto result = q_lx_val +
                                            lc * (vc_type(X[m].data() + i) + xloc_tmp) *
                                                vc_type(Q[sx_i + q][d].data() + i);
                                        result.store(Q[lx_i + n][d].data() + i);*/
                                        const auto result = q_lx_val +
                                            lc * (vc_type(X[m].data() + i) + xloc_tmp) *
                                                vc_type(combined_q + start_index_sx + q_i);
                                        result.store(combined_q + start_index_lx_n + q_i);
                                    }
                                }
                            }
                        }
                    }
                }
                for (int j = 0; j < geo.H_NX_XM4; j++) {
                    for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma GCC ivdep
                        for (int l = 0; l < geo.H_NX_ZM4; l++) {
                            const int q_i = to_q_index(j, k, l);
                            const auto rho = combined_q[start_index_rho + q_i];
                            combined_q[start_index_lx_n + q_i] *= rho;
                        }
                    }
                }
            }
            for (int dim = 0; dim < NDIM; dim++) {
								const int start_index_sx_d = (sx_i + dim) * q_face_offset + d * q_dir_offset;
                for (int j = 0; j < geo.H_NX_XM4; j++) {
                    for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                        for (int l = 0; l < geo.H_NX_ZM4; l++) {
                            const int q_i = to_q_index(j, k, l);
                            //const auto rho = Q[rho_i][d][i];
                            const auto rho = combined_q[start_index_rho + q_i];
                            //auto& v = Q[sx_i + dim][d][i];
                            auto& v = combined_q[start_index_sx_d + q_i];
                            combined_q[start_index_egas] += 0.5 * v * v * rho;
                            v *= rho;
                        }
                    }
                }
            }
            for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM4; l++) {
                        const int q_i = to_q_index(j, k, l);
                        const auto rho = combined_q[start_index_rho + q_i];
                        combined_q[start_index_pot + q_i] *= rho;
                        //Q[pot_i][d][i] *= rho;
                    }
                }
            }
            for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM4; l++) {
                        const int q_i = to_q_index(j, k, l);
                        //const auto rho = Q[rho_i][d][i];
                        const auto rho = combined_q[start_index_rho + q_i];
                        safe_real w = 0.0;
                        for (int si = 0; si < n_species_; si++) {
			                      const int start_index_sp_i = (spc_i + si) * q_face_offset + d * q_dir_offset;
                            w += combined_q[start_index_sp_i + q_i];
                            combined_q[start_index_sp_i + q_i] *= rho;
                            //w += Q[spc_i + si][d][i];
                            //Q[spc_i + si][d][i] *= rho;
                        }
                        w = 1.0 / w;
                        for (int si = 0; si < n_species_; si++) {
			                      const int start_index_sp_i = (spc_i + si) * q_face_offset + d * q_dir_offset;
                            combined_q[start_index_sp_i + q_i] *= w;
                            //Q[spc_i + si][d][i] *= w;
                        }
                    }
                }
            }
        }
    }
    return;
}

void post_recon_experimental(std::vector<std::vector<std::vector<safe_real>>>& Q,
    const hydro::x_type X, safe_real omega, bool angmom) {
    PROFILE();
    static const cell_geometry<NDIM, INX> geo;
    static const auto indices = geo.find_indices(2, geo.H_NX - 2);
    const auto dx = X[0][geo.H_DNX] - X[0][0];
    const auto xloc = geo.xloc();
    static constexpr auto levi_civita = geo.levi_civita();
    const int n_species_ = physics<NDIM>::get_n_species();
    auto dir = geo.direction();

    for (int d = 0; d < geo.NDIR; d++) {
        if (d != geo.NDIR / 2) {
            for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM4; l++) {
                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                        Q[sx_i][d][i] -= omega * (X[1][i] + 0.5 * xloc[d][1] * dx);
                        Q[sy_i][d][i] += omega * (X[0][i] + 0.5 * xloc[d][0] * dx);
                    }
                }
            }

            for (int n = 0; n < geo.NANGMOM; n++) {
                for (int q = 0; q < NDIM; q++) {
                    for (int m = 0; m < NDIM; m++) {
                        const auto lc = levi_civita[n][m][q];
                        if (lc != 0) {
                            const vc_type xloc_tmp = vc_type(0.5 * xloc[d][m] * dx);
                            for (int j = 0; j < geo.H_NX_XM4; j++) {
                                for (int k = 0; k < geo.H_NX_YM4; k++) {
                                    // TODO Reconstruct hotspot 2
                                    // for (int l = 0; l < geo.H_NX_ZM4; l++) {
                                    for (int l = 0; l < geo.H_NX_ZM4; l += vc_type::size()) {
                                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                                        const vc_type rho(Q[rho_i][d].data() + i);
                                        const vc_type q_lx_val(Q[lx_i + n][d].data() + i);
                                        // Q[lx_i + n][d][i] += lc *
                                        //   (X[m][i] + 0.5 * xloc[d][m] * dx) * Q[sx_i + q][d][i];
                                        const auto result = q_lx_val +
                                            lc * (vc_type(X[m].data() + i) + xloc_tmp) *
                                                vc_type(Q[sx_i + q][d].data() + i);
                                        result.store(Q[lx_i + n][d].data() + i);
                                    }
                                }
                            }
                        }
                    }
                }
                for (int j = 0; j < geo.H_NX_XM4; j++) {
                    for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma GCC ivdep
                        for (int l = 0; l < geo.H_NX_ZM4; l++) {
                            const int i = geo.to_index(j + 2, k + 2, l + 2);
                            const auto rho = Q[rho_i][d][i];
                            Q[lx_i + n][d][i] *= rho;
                        }
                    }
                }
            }
            for (int dim = 0; dim < NDIM; dim++) {
                for (int j = 0; j < geo.H_NX_XM4; j++) {
                    for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                        for (int l = 0; l < geo.H_NX_ZM4; l++) {
                            const int i = geo.to_index(j + 2, k + 2, l + 2);
                            const auto rho = Q[rho_i][d][i];
                            auto& v = Q[sx_i + dim][d][i];
                            Q[egas_i][d][i] += 0.5 * v * v * rho;
                            v *= rho;
                        }
                    }
                }
            }
            for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM4; l++) {
                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                        const auto rho = Q[rho_i][d][i];
                        Q[pot_i][d][i] *= rho;
                    }
                }
            }
            for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
#pragma ivdep
                    for (int l = 0; l < geo.H_NX_ZM4; l++) {
                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                        const auto rho = Q[rho_i][d][i];
                        safe_real w = 0.0;
                        for (int si = 0; si < n_species_; si++) {
                            w += Q[spc_i + si][d][i];
                            Q[spc_i + si][d][i] *= rho;
                        }
                        w = 1.0 / w;
                        for (int si = 0; si < n_species_; si++) {
                            Q[spc_i + si][d][i] *= w;
                        }
                    }
                }
            }
            /*for (int j = 0; j < geo.H_NX_XM4; j++) {
                for (int k = 0; k < geo.H_NX_YM4; k++) {
                    for (int l = 0; l < geo.H_NX_ZM4; l += vc_type::size()) {
                        const int i = geo.to_index(j + 2, k + 2, l + 2);
                        const vc_type rho(Q[rho_i][d].data() + i);
                        for (int n = 0; n < geo.NANGMOM; n++) {
                            const vc_type lx_i_tmp = vc_type(Q[lx_i + n][d].data() + i) * rho;
                            lxi_i_tmp.store(Q[lx_i + n][d].data() + i);
                            //Q[lx_i + n][d][i] *= rho;
                        }
                        for (int dim = 0; dim < NDIM; dim++) {
                            auto& v = Q[sx_i + dim][d][i];
                            Q[egas_i][d][i] += 0.5 * v * v * rho;
                            v *= rho;
                        }
                        Q[pot_i][d][i] *= rho;
                        safe_real w = 0.0;
                        for (int si = 0; si < n_species_; si++) {
                            w += Q[spc_i + si][d][i];
                            Q[spc_i + si][d][i] *= rho;
                        }
                        w = 1.0 / w;
                        for (int si = 0; si < n_species_; si++) {
                            Q[spc_i + si][d][i] *= w;
                        }
                    }
                }
            }*/
        }
    }
}
#pragma GCC pop_options
