/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "config.h"

/* Some standard headers. */
#include <fenv.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

#if defined(COOLING_EAGLE) && defined(CHEMISTRY_EAGLE) && defined(GADGET2_SPH)
#include "cooling/EAGLE/cooling_rates.h"
#include "cooling/EAGLE/cooling_tables.h"

/**
 * @brief Assign particle density and entropy corresponding to the
 * hydrogen number density and internal energy specified.
 *
 * @param p Particle data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 * @param internal_const Physical constants data structure
 * @param Z metallicity
 * @param nh Hydrogen number density (cgs units)
 * @param u Internal energy (cgs units)
 */
void set_quantities(struct part *restrict p, struct xpart *restrict xp,
                    const struct unit_system *restrict us,
                    const struct cooling_function_data *restrict cooling,
                    const struct cosmology *restrict cosmo,
                    const struct phys_const *restrict internal_const,
		    float Z, float nh, double u) {

  double hydrogen_number_density =
      nh * pow(units_cgs_conversion_factor(us, UNIT_CONV_LENGTH), 3);
  p->rho = hydrogen_number_density * internal_const->const_proton_mass /
           p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  float pressure = (u * cosmo->a * cosmo->a) *
                   cooling->internal_energy_from_cgs * p->rho *
                   (hydro_gamma_minus_one);
  p->entropy = pressure * (pow(p->rho, -hydro_gamma));
  xp->entropy_full = p->entropy;
  hydro_set_physical_internal_energy_dt(p,cosmo,0);

  p->chemistry_data.smoothed_metal_mass_fraction_total = Z;

  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    if (elem == 0) p->chemistry_data.smoothed_metal_mass_fraction[elem] = 0.752 - 3.29*Z;
    else if (elem == 1) p->chemistry_data.smoothed_metal_mass_fraction[elem] = 0.248 + 2.29*Z;
    else p->chemistry_data.smoothed_metal_mass_fraction[elem] *= Z;
  }

}

/**
 * @brief Produces contributions to cooling rates for different
 * hydrogen number densities, from different metals,
 * tests 1d and 4d table interpolations produce
 * same results for cooling rate, dlambda/du and temperature.
 */
int main(int argc, char **argv) {
  // Declare relevant structs
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const internal_const;
  struct hydro_props hydro_properties;
  struct entropy_floor_properties floor_props;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  struct space s;
  const char *parametersFileName = "./cooling_grid.yml";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  // Set some default values
  float redshift = 0.0, log_10_nh = -1;

  // Read options
  int param;
  while ((param = getopt(argc, argv, "z:d:")) != -1) switch (param) {
      case 'z':
        // read redshift
        redshift = atof(optarg);
        break;
      case 'd':
        // read log10 of hydrogen number density
        log_10_nh = atof(optarg);
        break;
      case '?':
        if (optopt == 'z')
          printf("Option -%c requires an argument.\n", optopt);
        else
          printf("Unknown option character `\\x%x'.\n", optopt);
        error("invalid option(s) to cooling_grid");
    }

  // Read the parameter file
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  // Init units
  units_init_from_params(&us, params, "InternalUnitSystem");

  // Init physical constants
  phys_const_init(&us, params, &internal_const);

  // Init porperties of hydro
  hydro_props_init(&hydro_properties, &internal_const, &us, params);

  // Init chemistry
  chemistry_init(params, &us, &internal_const, &chem_data);
  chemistry_first_init_part(&internal_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_print(&chem_data);

  // Init cosmology
  cosmology_init(params, &us, &internal_const, &cosmo);

  // Set redshift and associated quantities
  const float scale_factor = 1.0 / (1.0 + redshift);
  integertime_t ti_current =
      log(scale_factor / cosmo.a_begin) / cosmo.time_base;
  cosmology_update(&cosmo, &internal_const, ti_current);
  message("Redshift is %f", cosmo.z);

  // Init entropy floor
  entropy_floor_init(&floor_props,&internal_const,&us,&hydro_properties,params);

  // Init cooling
  cooling_init(params, &us, &internal_const, &hydro_properties, &cooling);
  cooling.H_reion_done = 1;
  cooling_print(&cooling);
  cooling_update(&cosmo, &cooling, &s);

  // Copy over the raw metals into the smoothed metals
  memcpy(&p.chemistry_data.smoothed_metal_mass_fraction,
         &p.chemistry_data.metal_mass_fraction,
         chemistry_element_count * sizeof(float));
  p.chemistry_data.smoothed_metal_mass_fraction_total =
      p.chemistry_data.metal_mass_fraction_total;

  // Calculate abundance ratios
  float abundance_ratio[(chemistry_element_count + 2)];
  abundance_ratio_to_solar(&p, &cooling, abundance_ratio);

  // open file
  FILE *output_iterations_file = fopen("cooling_iterations.dat", "w");
  if (output_iterations_file == NULL) {
    error("Error opening output iterations file!\n");
  }
  FILE *output_energy_file = fopen("cooling_energy_evolution.dat", "w");
  if (output_energy_file == NULL) {
    error("Error opening output energy file!\n");
  }

  // Set dts
  const float dt = 1.e-4;
  const float dt_therm = 1.e-4;

  // Define bounds for internal energy, density, metallicity
  const float log10_u_min_cgs = 10.f;
  const float log10_u_max_cgs = 18.f;
  const float log10_nh_min_cgs = -4.f;
  const float log10_nh_max_cgs = 0.f;
  const float log10_Z_min = -6.f;
  const float log10_Z_max = -2.f;
  const int n_u = 100;
  const int n_nh = 100;
  const int n_Z = 100;

  // Loop over internal energy
  for (int u_i = 0; u_i < n_u; u_i++) {
    double u_ini_cgs = exp10(log10_u_min_cgs + u_i * (log10_u_max_cgs - log10_u_min_cgs)/n_u);
    double u_ini = u_ini_cgs/units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY_PER_UNIT_MASS);

    // Loop over hydrogen number density
    for (int nh_i = 0; nh_i < n_nh; nh_i++) {
      float nh_cgs = exp10(log10_nh_min_cgs + nh_i * (log10_nh_max_cgs - log10_nh_min_cgs)/n_nh);

      // Loop over metallicities
      for (int Z_i = 0; Z_i < n_Z; Z_i++) {
        float Z = exp10(log10_Z_min + Z_i * (log10_Z_max - log10_Z_min)/n_Z);
  
	// Update particle data
        set_quantities(&p, &xp, &us, &cooling, &cosmo, &internal_const, Z,
		       nh_cgs, u_ini_cgs);

        // Cool the particle
	cooling_cool_part(&internal_const, &us, &cosmo,
                       &hydro_properties, &floor_props,
                       // Set to non const for counting, remove for production
                       //const struct cooling_function_data *cooling,
                       &cooling, &p, &xp, dt, dt_therm);
        
	// Update the particle's internal energy
	const float du_dt_new = hydro_get_physical_internal_energy_dt(&p, &cosmo);
	hydro_set_physical_internal_energy(&p, &xp, &cosmo, u_ini + dt * du_dt_new);

	// Get the final energy of the particle
	const double u_final_cgs = hydro_get_physical_internal_energy(&p,&xp,&cosmo) *
				   units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY_PER_UNIT_MASS);

	// Save relevant data
	fprintf(output_iterations_file, "%.5e %d %d %d\n", 
		u_ini_cgs, cooling.bisection_cooling_bound_iterations, 
		cooling.bisection_heating_bound_iterations, 
		cooling.bisection_iterations);
	
	fprintf(output_energy_file, "%.5e %.5e\n", u_ini_cgs, u_final_cgs);

	// Zero counters
	cooling.bisection_cooling_bound_iterations = 0;
	cooling.bisection_heating_bound_iterations = 0;
	cooling.bisection_iterations = 0;
        
      }
    }
  }

  fclose(output_iterations_file);
  fclose(output_energy_file);
  message("done cooling rates test");

  /* Clean everything */
  cosmology_clean(&cosmo);
  cooling_clean(&cooling);

  free(params);
  return 0;
}

#else

int main(int argc, char **argv) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  message("This test is only defined for the EAGLE cooling model.");
  return 0;
}
#endif
