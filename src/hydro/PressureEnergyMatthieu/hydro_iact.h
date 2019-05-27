/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_MINIMAL_HYDRO_IACT_H
#define SWIFT_MINIMAL_HYDRO_IACT_H

/**
 * @file Minimal/hydro_iact.h
 * @brief Minimal conservative implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term without switches is implemented. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adiabatic_index.h"
#include "minmax.h"

#include "./hydro_parameters.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Get the internal energies */
  const float ui = pi->u;
  const float uj = pj->u;
  
  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float r_over_hi = r * hi_inv;
  kernel_deval(r_over_hi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + r_over_hi * wi_dx);

  pi->rho += mj * wi;
  
  pi->P_bar += mj * uj * wi;
  pi->P_bar_dh -= mj * uj * (hydro_dimension * wi + r_over_hi * wi_dx);
  
  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float r_over_hj = r * hj_inv;
  kernel_deval(r_over_hj, &wj, &wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + r_over_hj * wj_dx);

  pj->rho += mi * wj;
  
  pj->P_bar += mi * ui * wj;
  pj->P_bar_dh -= mi * ui * (hydro_dimension * wj + r_over_hj * wj_dx);
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

  float wi, wi_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get the internal energies */
  const float uj = pj->u;
  
  /* Get r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  const float h_inv = 1.f / hi;
  const float r_over_hi = r * h_inv;
  kernel_deval(r_over_hi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + r_over_hi * wi_dx);

  pi->rho += mj * wi;
  
  pi->P_bar += mj * uj * wi;
  pi->P_bar_dh -= mj * uj * (hydro_dimension * wi + r_over_hi * wi_dx);
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and r inverse. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float r_over_hi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(r_over_hi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float r_over_hj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(r_over_hj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;
  
  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float ui = pi->u;
  const float uj = pj->u;
  const float P_bar_i = pi->P_bar;
  const float P_bar_j = pj->P_bar;
  const float P_bar_dh_i = pi->P_bar_dh;
  const float P_bar_dh_j = pj->P_bar_dh;
  const float n_i = pi->density.wcount; // Note this is legal since we don't have a union
  const float n_j = pj->density.wcount; // Note this is legal since we don't have a union
  const float n_dh_i = pi->density.wcount_dh;
  const float n_dh_j = pj->density.wcount_dh;
  
  /* Compute f_ij terms */
  const float bracket1_i = hi * P_bar_dh_i / (hydro_dimension * n_i);
  const float bracket1_j = hj * P_bar_dh_j / (hydro_dimension * n_j);

  const float bracket2_i = 1.f + (hi *  n_dh_i / (hydro_dimension * n_i));
  const float bracket2_j = 1.f + (hj *  n_dh_j / (hydro_dimension * n_j));

  const float f_ij = 1.f - (1.f / (hydro_gamma_minus_one * mj * uj) ) * bracket1_i / bracket2_i;
  const float f_ji = 1.f - (1.f / (hydro_gamma_minus_one * mi * ui) ) * bracket1_j / bracket2_j;
  
  /* SPH acceleration term (we multiply by the mass at the end) */
  const float sph_acc_term =
    hydro_gamma_minus_one * hydro_gamma_minus_one * ui * uj *
    ((f_ij / P_bar_i) * wi_dr + (f_ji / P_bar_j) * wj_dr) * r_inv;
  
  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = sqrtf(ui * hydro_gamma_minus_one * hydro_gamma);
  const float cj = sqrtf(uj * hydro_gamma_minus_one * hydro_gamma);
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;
  
  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = hydro_gamma_minus_one * hydro_gamma_minus_one * ui * uj * (f_ij / P_bar_i) * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = hydro_gamma_minus_one * hydro_gamma_minus_one * ui * uj * (f_ji / P_bar_j) * dvdr * r_inv * wj_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {


    /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and r inverse. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float r_over_hi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(r_over_hi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float r_over_hj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(r_over_hj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;
  
  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float ui = pi->u;
  const float uj = pj->u;
  const float P_bar_i = pi->P_bar;
  const float P_bar_j = pj->P_bar;
  const float P_bar_dh_i = pi->P_bar_dh;
  const float P_bar_dh_j = pj->P_bar_dh;
  const float n_i = pi->density.wcount; // Note this is legal since we don't have a union
  const float n_j = pj->density.wcount; // Note this is legal since we don't have a union
  const float n_dh_i = pi->density.wcount_dh;
  const float n_dh_j = pj->density.wcount_dh;
  
  /* Compute f_ij terms */
  const float bracket1_i = hi * P_bar_dh_i / (hydro_dimension * n_i);
  const float bracket1_j = hj * P_bar_dh_j / (hydro_dimension * n_j);

  const float bracket2_i = 1.f + (hi *  n_dh_i / (hydro_dimension * n_i));
  const float bracket2_j = 1.f + (hj *  n_dh_j / (hydro_dimension * n_j));

  const float f_ij = 1.f - (1.f / (hydro_gamma_minus_one * mj * uj) ) * bracket1_i / bracket2_i;
  const float f_ji = 1.f - (1.f / (hydro_gamma_minus_one * mi * ui) ) * bracket1_j / bracket2_j;
  
  /* SPH acceleration term (we multiply by the mass at the end) */
  const float sph_acc_term =
    hydro_gamma_minus_one * hydro_gamma_minus_one * ui * uj *
    ((f_ij / P_bar_i) * wi_dr + (f_ji / P_bar_j) * wj_dr) * r_inv;
  
  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = sqrtf(ui * hydro_gamma_minus_one * hydro_gamma);
  const float cj = sqrtf(uj * hydro_gamma_minus_one * hydro_gamma);
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;
  
  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = hydro_gamma_minus_one * hydro_gamma_minus_one * ui * uj * (f_ij / P_bar_i) * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
}

/**
 * @brief Timestep limiter loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_limiter(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  /* Nothing to do here if both particles are active */
}

/**
 * @brief Timestep limiter loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_limiter(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  /* Wake up the neighbour? */
  if (pi->force.v_sig > const_limiter_max_v_sig_ratio * pj->force.v_sig) {

    pj->wakeup = time_bin_awake;
  }
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
