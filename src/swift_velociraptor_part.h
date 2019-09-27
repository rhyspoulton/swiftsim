/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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
#ifndef SWIFT_VELOCIRAPTOR_PART_H
#define SWIFT_VELOCIRAPTOR_PART_H

#include "part_type.h"
#include "chemistry.h"

/**
 * @brief SWIFT/VELOCIraptor particle.
 *
 * This should match the structure Swift::swift_vel_part
 * defined in the file NBodylib/src/NBody/SwiftParticle.h
 * of the VELOCIraptor code.
 */
struct swift_vel_part {

  /*! Particle ID. */
  long long id;

  /*! Particle position. */
  double x[3];

  /*! Particle velocity. */
  float v[3];

  #ifndef VR_NOMASS
  /*! Particle mass. */
  float mass;
  #endif

  /*! Gravitational potential */
  float potential;

  /*! Internal energy of gas particle */
  float u;

  /*! Temperature of a gas particle */
  float T;

  /*! Type of the #gpart (DM, gas, star, ...) */
  enum part_type type;

  /*! MPI rank on which this #gpart lives on the SWIFT side. */
  int task;

  /*! Index of this #gpart in the global array of this rank on the SWIFT
    side. */
  int index;
};

/* SWIFT/VELOCIraptor chemistry data. */
struct swift_vel_chemistry_data {

    /*! Fraction of the particle mass in a given element */
    float metal_mass_fraction[chemistry_element_count];

    /*! Fraction of the particle mass in *all* metals */
    float metal_mass_fraction_total;

    /*! Mass coming from SNIa */
    float mass_from_SNIa;

    /*! Fraction of total gas mass in metals coming from SNIa */
    float metal_mass_fraction_from_SNIa;

    /*! Mass coming from AGB */
    float mass_from_AGB;

    /*! Fraction of total gas mass in metals coming from AGB */
    float metal_mass_fraction_from_AGB;

    /*! Mass coming from SNII */
    float mass_from_SNII;

    /*! Fraction of total gas mass in metals coming from SNII */
    float metal_mass_fraction_from_SNII;

    /*! Fraction of total gas mass in Iron coming from SNIa */
    float iron_mass_fraction_from_SNIa;

};

/* SWIFT/VELOCIraptor chemistry of black holes. */
struct swift_vel_chemistry_bh_data {

    /*! Mass in a given element */
    float metal_mass[chemistry_element_count];

    /*! Mass in *all* metals */
    float metal_mass_total;

    /*! Mass coming from SNIa */
    float mass_from_SNIa;

    /*! Mass coming from AGB */
    float mass_from_AGB;

    /*! Mass coming from SNII */
    float mass_from_SNII;

    /*! Metal mass coming from SNIa */
    float metal_mass_from_SNIa;

    /*! Metal mass coming from AGB */
    float metal_mass_from_AGB;

    /*! Metal mass coming from SNII */
    float metal_mass_from_SNII;

    /*! Iron mass coming from SNIa */
    float iron_mass_from_SNIa;
};


/* SWIFT/VELOCIraptor gas particle. */
struct swift_vel_gas_part {
    /*! Particle smoothing length. */
    float h;

    /*! Particle internal energy. */
    float u;

    /*! Time derivative of the internal energy. */
    float u_dt;

    /*! Particle density. */
    float rho;

    /*! Particle pressure (weighted) */
    float pressure_bar;

    /* Store viscosity information in a separate struct. */
    struct {

    /*! Particle velocity divergence */
    float div_v;

    /*! Particle velocity divergence from previous step */
    float div_v_previous_step;

    /*! Artificial viscosity parameter */
    float alpha;

    /*! Signal velocity */
    float v_sig;

    } viscosity;

    /* Store thermal diffusion information in a separate struct. */
    struct {

    /*! del^2 u, a smoothed quantity */
    float laplace_u;

    /*! Thermal diffusion coefficient */
    float alpha;

    } diffusion;

    /* Chemistry information */
    struct swift_vel_chemistry_data chemistry_data;

    /*! swift index */
    int index;
};


/* SWIFT/VELOCIraptor star particle. */
struct swift_vel_star_part {

    /*! Birth time (or scalefactor)*/
    float birth_time;

    /*! Birth density */
    float birth_density;

    /*! Birth temperature */
    float birth_temperature;

    /*! Feedback energy fraction */
    float f_E;

    /*! Chemistry structure */
    struct swift_vel_chemistry_data chemistry_data;

    /*! swift index */
    int index;
};

/* SWIFT/VELOCIraptor black hole particle. */
struct swift_vel_bh_part {

    /*! Formation time (or scale factor)*/
    float formation_time;

    /*! Subgrid mass of the black hole */
    float subgrid_mass;

    /*! Total accreted mass of the black hole (including accreted mass onto BHs
    * that were merged) */
    float total_accreted_mass;

    /*! Energy reservoir for feedback */
    float energy_reservoir;

    /*! Instantaneous accretion rate */
    float accretion_rate;

    /*! Density of the gas surrounding the black hole. */
    float rho_gas;

    /*! Smoothed sound speed of the gas surrounding the black hole. */
    float sound_speed_gas;

    /*! Smoothed velocity (peculiar) of the gas surrounding the black hole */
    float velocity_gas[3];

    /*! Curl of the velocity field around the black hole */
    float circular_velocity_gas[3];

    /*! Number of seeds in this BH (i.e. itself + the merged ones) */
    int cumulative_number_seeds;

    /*! Total number of BH merger events (i.e. not including all progenies) */
    int number_of_mergers;

    /*! Chemistry information (e.g. metal content at birth, swallowed metal
    * content, etc.) */
    struct swift_vel_chemistry_bh_data chemistry_data;

    /*! swift index */
    int index;
};

#endif /* SWIFT_VELOCIRAPTOR_PART_H */
