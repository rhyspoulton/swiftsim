/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_VELOCIRAPTOR_IO_H
#define SWIFT_VELOCIRAPTOR_IO_H

/* Config parameters. */
#include "../config.h"
#include "velociraptor_struct.h"

INLINE static void velociraptor_convert_part_groupID(const struct engine* e,
                                                     const struct part* p,
                                                     const struct xpart* xp,
                                                     long long* ret) {
  if (p->gpart == NULL)
    ret[0] = 0.f;
  else {
    const ptrdiff_t offset = p->gpart - e->s->gparts;
    *ret = (e->s->gpart_group_data + offset)->groupID;
  }
}

INLINE static void velociraptor_convert_spart_groupID(const struct engine* e,
                                                      const struct spart* sp,
                                                      long long* ret) {
  if (sp->gpart == NULL)
    ret[0] = 0.f;
  else {
    const ptrdiff_t offset = sp->gpart - e->s->gparts;
    *ret = (e->s->gpart_group_data + offset)->groupID;
  }
}

INLINE static void velociraptor_convert_bpart_groupID(const struct engine* e,
                                                      const struct bpart* bp,
                                                      long long* ret) {
  if (bp->gpart == NULL)
    ret[0] = 0.f;
  else {
    const ptrdiff_t offset = bp->gpart - e->s->gparts;
    *ret = (e->s->gpart_group_data + offset)->groupID;
  }
}

__attribute__((always_inline)) INLINE static int velociraptor_write_parts(
    const struct part* parts, const struct xpart* xparts,
    struct io_props* list) {

  list[0] = io_make_output_field_convert_part(
      "GroupID", LONGLONG, 1, UNIT_CONV_NO_UNITS, parts, xparts,
      velociraptor_convert_part_groupID);

  return 1;
}

__attribute__((always_inline)) INLINE static int velociraptor_write_gparts(
    const struct velociraptor_gpart_data* group_data, struct io_props* list) {

  list[0] = io_make_output_field("GroupID", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 group_data, groupID);

  return 1;
}

__attribute__((always_inline)) INLINE static int velociraptor_write_sparts(
    const struct spart* sparts, struct io_props* list) {

  list[0] = io_make_output_field_convert_spart(
      "GroupID", LONGLONG, 1, UNIT_CONV_NO_UNITS, sparts,
      velociraptor_convert_spart_groupID);

  return 1;
}

__attribute__((always_inline)) INLINE static int velociraptor_write_bparts(
    const struct bpart* bparts, struct io_props* list) {

  list[0] = io_make_output_field_convert_bpart(
      "GroupID", LONGLONG, 1, UNIT_CONV_NO_UNITS, bparts,
      velociraptor_convert_bpart_groupID);

  return 1;
}


/**
 * @brief Copy every velociraptor defined most bound particle #part into the parts_written array.
 *
 * @param ibmp whether select only most bound particles or any particle in a group
 * @param groupinfo group information used to select particles
 * @param parts The array of #part containing all particles.
 * @param xparts The array of #xpart containing all particles.
 * @param parts_written The array of #part to fill with particles we want to
 * write.
 * @param xparts_written The array of #xpart  to fill with particles we want to
 * write.
 * @param Nparts The total number of #part.
 * @param Nparts_written The total number of #part to write.
 */
void io_collect_velociraptor_parts_to_write(const int imbp,
                                const struct velociraptor_gpart_data* restrict groupinfo,
                                const struct part* restrict parts,
                                const struct xpart* restrict xparts,
                                struct part* restrict parts_written,
                                struct xpart* restrict xparts_written,
                                const size_t Nparts,
                                const size_t Nparts_written) {

  size_t count = 0;
  /* Loop over all parts */
  if (imbp) {
      for (size_t i = 0; i < Nparts; ++i) {
        //if particle is every the most bound particle of a group, write it.
        if (groupinfo[i].groupID_lastmbp>0) {
          parts_written[count] = parts[i];
          xparts_written[count] = xparts[i];
          count++;
        }
      }
  }
  else {
      for (size_t i = 0; i < Nparts; ++i) {

        if (groupinfo[i].groupID>0) {

          parts_written[count] = parts[i];
          xparts_written[count] = xparts[i];
          count++;
        }
      }
  }

  /* Check that everything is fine */
  if (count != Nparts_written)
    error("Collected the wrong number of particles (%zu vs. %zu expected)",
          count, Nparts_written);
}

/**
 * @brief Copy every velociraptor defined most bound particle #spart into the sparts_written array.
 *
 * @param ibmp whether select only most bound particles or any particle in a group
 * @param groupinfo group information used to select particles
 * @param sparts The array of #spart containing all particles.
 * @param sparts_written The array of #spart to fill with particles we want to
 * write.
 * @param Nsparts The total number of #spart.
 * @param Nsparts_written The total number of #spart to write.
 */

void io_collect_velociraptor_sparts_to_write(const int imbp,
                                const struct velociraptor_gpart_data* restrict groupinfo,
                                const struct spart* restrict sparts,
                                struct spart* restrict sparts_written,
                                const size_t Nsparts,
                                const size_t Nsparts_written) {

  size_t count = 0;

  /* Loop over all parts */
  if (imbp) {
      for (size_t i = 0; i < Nsparts; ++i) {
          //if particle is every the most bound particle of a group, write it.
          if (groupinfo[i].groupID_lastmbp>0) {
          sparts_written[count] = sparts[i];
          count++;
        }
      }
  }
  else {
      for (size_t i = 0; i < Nsparts; ++i) {

        if (groupinfo[i].groupID>0) {
            sparts_written[count] = sparts[i];
            count++;
        }
      }
  }

  /* Check that everything is fine */
  if (count != Nsparts_written)
    error("Collected the wrong number of s-particles (%zu vs. %zu expected)",
          count, Nsparts_written);
}

/**
 * @brief Copy every velociraptor defined most bound particle #bpart into the bparts_written array.
 *
 * @param ibmp whether select only most bound particles or any particle in a group
 * @param groupinfo group information used to select particles
 * @param bparts The array of #bpart containing all particles.
 * @param bparts_written The array of #bpart to fill with particles we want to
 * write.
 * @param Nbparts The total number of #bpart.
 * @param Nbparts_written The total number of #bpart to write.
 */

void io_collect_velociraptor_bparts_to_write(const int imbp,
                                const struct velociraptor_gpart_data* restrict groupinfo,
                                const struct bpart* restrict bparts,
                                struct bpart* restrict bparts_written,
                                const size_t Nbparts,
                                const size_t Nbparts_written) {

  size_t count = 0;

  /* Loop over all parts */
  if (imbp) {
      for (size_t i = 0; i < Nbparts; ++i) {
          //if particle is every the most bound particle of a group, write it.
          if (groupinfo[i].groupID_lastmbp>0) {
          bparts_written[count] = bparts[i];
          count++;
        }
      }
  }
  else {
      for (size_t i = 0; i < Nbparts; ++i) {
        if (groupinfo[i].groupID>0) {
            bparts_written[count] = bparts[i];
            count++;
        }
      }
  }

  /* Check that everything is fine */
  if (count != Nbparts_written)
    error("Collected the wrong number of b-particles (%zu vs. %zu expected)",
          count, Nbparts_written);
}

/**
 * @brief Copy every velociraptor defined most bound particle #gpart into the gparts_written array.
 *
 * @param ibmp whether select only most bound particles or any particle in a group
 * @param groupinfo group information used to select particles
 * @param gparts The array of #gpart containing all particles.
 * @param gparts_written The array of #gpart to fill with particles we want to
 * write.
 * @param Ngparts The total number of #gpart.
 * @param Ngparts_written The total number of #gpart to write.
 */
void io_collect_velociraptor_gparts_to_write(const int imbp,
    const struct velociraptor_gpart_data* restrict groupinfo,
    const struct gpart* restrict gparts,
    const struct velociraptor_gpart_data* restrict vr_data,
    struct gpart* restrict gparts_written,
    const size_t Ngparts, const size_t Ngparts_written) {

  size_t count = 0;

  /* Loop over all parts */
  if (imbp) {
      for (size_t i = 0; i < Ngparts; ++i) {
          //if particle is every the most bound particle of a group, write it.
          if (groupinfo[i].groupID_lastmbp>0) {
          gparts_written[count] = gparts[i];
          count++;
        }
    }
  }
  else {
      for (size_t i = 0; i < Ngparts; ++i) {

        if (groupinfo[i].groupID>0) {
            gparts_written[count] = gparts[i];
            count++;
        }
      }
  }

  /* Check that everything is fine */
  if (count != Ngparts_written)
    error("Collected the wrong number of g-particles (%zu vs. %zu expected)",
          count, Ngparts_written);
}

#if defined(HAVE_HDF5) && defined(WITH_MPI) && defined(HAVE_PARALLEL_HDF5)

/**
 * @brief Prepares a velociraptor based snippet file for a parallel write
 *
 * @param e The #engine.
 * @param baseName The base name of the snapshots.
 * @param N_total The total number of particles of each type to write.
 * @param whether file is most bound particle set only.
 * @param internal_units The #unit_system used internally.
 * @param snapshot_units The #unit_system used in the snapshots.
 */
void prepare_velociraptor_snippet_file(struct engine* e, const char* baseName,
                  long long N_total[6],
                  const int ifile_mbp,
                  const struct unit_system* internal_units,
                  const struct unit_system* snapshot_units) {

  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;
  const struct gpart* gparts = e->s->gparts;
  const struct spart* sparts = e->s->sparts;
  const struct bpart* bparts = e->s->bparts;
  const struct groupinfo *group_info = e->s->gpart_group_data;
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Nblackholes = e->s->nr_bparts;
  const size_t Ndm = Ntot-Ngas-Nstars-Nblackholes;

  struct part* parts_written=NULL;
  struct xpart* xparts_written=NULL;
  struct gpart* gparts_written=NULL;
  struct spart* sparts_written=NULL;
  struct bpart* bparts_written=NULL;
  ptrdiff_t offset;
  int counter[swift_type_count] = {0};

  struct swift_params* params = e->parameter_file;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;

  FILE* xmfFile = 0;
  int numFiles = 1;

  /* First time, we need to create the XMF file */
  if (e->snapshot_output_count == 0) xmf_create_file(baseName);

  /* Prepare the XMF file for the new entry */
  xmfFile = xmf_prepare_file(baseName);

  /* HDF5 File name */
  char fileName[FILENAME_BUFFER_SIZE];
  if (e->snapshot_int_time_label_on)
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%06i.hdf5", baseName,
             (int)round(e->time));
  else
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%04i.hdf5", baseName,
             e->snapshot_output_count);

  /* Open HDF5 file with the chosen parameters */
  hid_t h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", fileName);

  /* Write the part of the XMF file corresponding to this
   * specific output */
  xmf_write_outputheader(xmfFile, fileName, e->time);

  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  hid_t h_grp =
      H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_TIME);
  const double factor_length =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  /* Print the relevant information and print status */
  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute(h_grp, "Time", DOUBLE, &dblTime, 1);
  const int dimension = (int)hydro_dimension;
  io_write_attribute(h_grp, "Dimension", INT, &dimension, 1);
  io_write_attribute(h_grp, "Redshift", DOUBLE, &e->cosmology->z, 1);
  io_write_attribute(h_grp, "Scale-factor", DOUBLE, &e->cosmology->a, 1);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  if (ifile_mbp) io_write_attribute_s(h_grp, "SnapshotFormat", "VR_MostBoundParticle_Snippet");
  else io_write_attribute_s(h_grp, "SnapshotFormat", "VR_InGroup_Snippet");
  time_t tm = time(NULL);
  io_write_attribute_s(h_grp, "Snapshot date", ctime(&tm));
  io_write_attribute_s(h_grp, "RunName", e->run_name);

  /* GADGET-2 legacy values */
  /* Number of particles of each type */
  long long N_Total_written[swift_type_count] = {0};
  unsigned int numParticles[swift_type_count] = {0};
  unsigned int numParticlesHighWord[swift_type_count] = {0};

  if (ifile_mbp) {
      for (int i=0; i<e->s->nr_gparts;i++) {
        if (group_info[i].groupID_lastmbp==0) continue;
        N_Total_written[gparts[i].type]++;
      }
  }
  else {
      for (int i=0; i<e->s->nr_gparts;i++) {
        if (group_info[i].groupID==0) continue;
        N_Total_written[gparts[i].type]++;
      }
  }
  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);
  }
  io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, N_total,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total", UINT, numParticles,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total_HighWord", UINT,
                     numParticlesHighWord, swift_type_count);
  double MassTable[6] = {0., 0., 0., 0., 0., 0.};
  io_write_attribute(h_grp, "MassTable", DOUBLE, MassTable, swift_type_count);
  unsigned int flagEntropy[swift_type_count] = {0};
  flagEntropy[0] = writeEntropyFlag();
  io_write_attribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                     swift_type_count);
  io_write_attribute(h_grp, "NumFilesPerSnapshot", INT, &numFiles, 1);

  /* Close header */
  H5Gclose(h_grp);

  /* Print the code version */
  io_write_code_description(h_file);

  /* Print the run's policy */
  io_write_engine_policy(h_file, e);

  /* Print the SPH parameters */
  if (e->policy & engine_policy_hydro) {
    h_grp = H5Gcreate(h_file, "/HydroScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating SPH group");
    hydro_props_print_snapshot(h_grp, e->hydro_properties);
    hydro_write_flavour(h_grp);
    H5Gclose(h_grp);
  }

  /* Print the subgrid parameters */
  h_grp = H5Gcreate(h_file, "/SubgridScheme", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating subgrid group");
  entropy_floor_write_flavour(h_grp);
  cooling_write_flavour(h_grp, e->cooling_func);
  chemistry_write_flavour(h_grp);
  tracers_write_flavour(h_grp);
  H5Gclose(h_grp);

  /* Print the gravity parameters */
  if (e->policy & engine_policy_self_gravity) {
    h_grp = H5Gcreate(h_file, "/GravityScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating gravity group");
    gravity_props_print_snapshot(h_grp, e->gravity_properties);
    H5Gclose(h_grp);
  }

  /* Print the stellar parameters */
  if (e->policy & engine_policy_stars) {
    h_grp = H5Gcreate(h_file, "/StarsScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating stars group");
    stars_props_print_snapshot(h_grp, e->stars_properties);
    H5Gclose(h_grp);
  }

  /* Print the cosmological parameters */
  h_grp =
      H5Gcreate(h_file, "/Cosmology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cosmology group");
  if (e->policy & engine_policy_cosmology)
    io_write_attribute_i(h_grp, "Cosmological run", 1);
  else
    io_write_attribute_i(h_grp, "Cosmological run", 0);
  cosmology_write_model(h_grp, e->cosmology);
  H5Gclose(h_grp);

  /* Print the runtime parameters */
  h_grp =
      H5Gcreate(h_file, "/Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, 1);
  H5Gclose(h_grp);

  /* Print the runtime unused parameters */
  h_grp = H5Gcreate(h_file, "/UnusedParameters", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, 0);
  H5Gclose(h_grp);

  /* Print the system of Units used in the spashot */
  io_write_unit_system(h_file, snapshot_units, "Units");

  /* Print the system of Units used internally */
  io_write_unit_system(h_file, internal_units, "InternalCodeUnits");


  io_collect_velociraptor_sparts_to_write(ifile_mbp, groupinfo, sparts, sparts_written, Nsparts, Nsparts_written);
  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (N_total[ptype] == 0) continue;
    if (N_total_written[ptype] == 0) continue;
    //allcoate mem to store information
    switch (ptype) {
        case swift_type_gas:
            offset = parts->gpart - e->s->gparts;
            parts_written=(struct part *)malloc(sizeof(struct part)*N_total_written[ptype]));
            xparts_written=(struct xpart *)malloc(sizeof(struct xpart)*N_total_written[ptype]));
            io_collect_velociraptor_parts_to_write(ifile_mbp, group_info+offset, parts, xparts, parts_written, xparts_written, Ngas, N_total_written[ptype]);
            break;
        case swift_type_dark_matter:
            offset = 0;
            gparts_written=(struct gpart *)malloc(sizeof(struct gpart)*N_total_written[ptype]));
            io_collect_velociraptor_gparts_to_write(ifile_mbp, group_info, gparts, gparts_written, Ndm, N_total_written[ptype]);
            break;
        case swift_type_stars:
            offset = sparts->gpart - e->s->gparts;
            sparts_written=(struct spart *)malloc(sizeof(struct spart)*N_total_written[ptype]));
            io_collect_velociraptor_sparts_to_write(ifile_mbp, group_info + offset, sparts, sparts_written, Nstar, N_total_written[ptype]);
            break;
        case swift_type_black_hole:
            offset = bparts->gpart - e->s->gparts;
            bparts_written=(struct bpart *)malloc(sizeof(struct bpart)*N_total_written[ptype]));
            io_collect_velociraptor_bparts_to_write(ifile_mbp, group_info + offset, bparts, bparts_written, Nblackholes, N_total_written[ptype]);
            break;
    }

    /* Add the global information for that particle type to
     * the XMF meta-file */
    xmf_write_groupheader(xmfFile, fileName, N_total[ptype],
                          (enum part_type)ptype);

    /* Create the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0)
      error("Error while opening particle group %s.", partTypeGroupName);

    int num_fields = 0;
    struct io_props list[100];

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas:
        hydro_write_particles(parts_written, xparts_written, list, &num_fields);
        num_fields += chemistry_write_particles(parts_written, list + num_fields);
        if (with_cooling || with_temperature) {
          num_fields += cooling_write_particles(
              parts_writen, xparts_written, list + num_fields, e->cooling_func);
        }
        num_fields += tracers_write_particles(parts_written, xparts_written, list + num_fields,
                                              with_cosmology);
        num_fields += star_formation_write_particles(parts_written, xparts_written, list + num_fields);
        num_fields += velociraptor_write_parts(parts_written, xparts_written, list + num_fields);
        break;

      case swift_type_dark_matter:
        darkmatter_write_particles(gparts_written, list, &num_fields);
        num_fields += velociraptor_write_gparts(gparts_written, list + num_fields);
        break;

      case swift_type_stars:
        stars_write_particles(sparts_written, list, &num_fields);
        num_fields += chemistry_write_sparticles(sparts_written, list + num_fields);
        num_fields += tracers_write_sparticles(sparts_written, list + num_fields, with_cosmology);
        break;

      case swift_type_black_hole:
        black_holes_write_particles(bparts_written, list, &num_fields);
        num_fields += chemistry_write_bparticles(bparts_written, list + num_fields);
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    free(parts_written);
    free(xparts_written);
    free(gparts_written);
    free(sparts_written);
    free(bparts_written);
    parts_written=NULL;
    xparts_written=NULL;
    gparts_written=NULL;
    xparts_written=NULL;
    bparts_written=NULL;

    /* Prepare everything that is not cancelled */
    for (int i = 0; i < num_fields; ++i) {

      /* Did the user cancel this field? */
      char field[PARSER_MAX_LINE_SIZE];
      sprintf(field, "SelectOutput:%.*s_%s", FIELD_BUFFER_SIZE, list[i].name,
              part_type_names[ptype]);
      int should_write = parser_get_opt_param_int(params, field, 1);

      if (should_write)
        prepareArray(e, h_grp, fileName, xmfFile, partTypeGroupName, list[i],
                     N_total[ptype], snapshot_units);
    }

    /* Close particle group */
    H5Gclose(h_grp);

    /* Close this particle group in the XMF file as well */
    xmf_write_groupfooter(xmfFile, (enum part_type)ptype);
  }

  /* Write LXMF file descriptor */
  xmf_write_outputfooter(xmfFile, e->snapshot_output_count, e->time);

  /* Close the file for now */
  H5Fclose(h_file);
}

/**
 * @brief Writes an parallel HDF5 snippet output file (GADGET-3 type) with
 *its XMF descriptor based on VR results.
 *
 * @param e The engine containing all the system.
 * @param baseName The common part of the snapshot file name.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 * @param mpi_rank The MPI rank of this node.
 * @param mpi_size The number of MPI ranks.
 * @param comm The MPI communicator.
 * @param info The MPI information object
 *
 * Creates an HDF5 output file and writes the particles
 * contained in the engine. If such a file already exists, it is
 * erased and replaced by the new one.
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_velociraptor_snippet_output_parallel(struct engine* e, const char* baseName,
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units,
                           int mpi_rank, int mpi_size, MPI_Comm comm,
                           MPI_Info info) {

  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;
  const struct gpart* gparts = e->s->gparts;
  const struct spart* sparts = e->s->sparts;
  const struct bpart* bparts = e->s->bparts;
  struct swift_params* params = e->parameter_file;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
  const int ifile_mbp = e->stf_output_mostboundparticles;
  ptrdiff_t offset;
  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Nblackholes = e->s->nr_bparts;
  // const size_t Nbaryons = Ngas + Nstars;
  // const size_t Ndm = Ntot > 0 ? Ntot - Nbaryons : 0;

  /* Number of particles that we will write */
  size_t Ntot_written, Ngas_written, Nstars_written, Nblackholes_written;
  if (ifile_mbp) {
      for (int i=0; i<e->s->nr_gparts;i++) {
        if (group_info[i].groupID_lastmbp==0) continue;
        if (gparts[i].type==swift_type_gas) Ngas_written++;
        else if (gparts[i].type==swift_type_stars) Nstars_written++;
        else if (gparts[i].type==swift_type_black_hole) Nblackholes_written++;
        Ntot_written++;

      }
  }
  else {
      for (int i=0; i<e->s->nr_gparts;i++) {
        if (group_info[i].groupID==0) continue;
        if (gparts[i].type==swift_type_gas) Ngas_written++;
        else if (gparts[i].type==swift_type_stars) Nstars_written++;
        else if (gparts[i].type==swift_type_black_hole) Nblackholes_written++;
        Ntot_written++;
      }
  }
  const size_t Nbaryons_written =
      Ngas_written + Nstars_written + Nblackholes_written;
  const size_t Ndm_written =
      Ntot_written > 0 ? Ntot_written - Nbaryons_written : 0;

  /* Compute offset in the file and total number of particles */
  size_t N[swift_type_count] = {Ngas_written,   Ndm_written,        0, 0,
                                Nstars_written, Nblackholes_written};
  long long N_total[swift_type_count] = {0};
  long long offset[swift_type_count] = {0};
  MPI_Exscan(&N, &offset, swift_type_count, MPI_LONG_LONG_INT, MPI_SUM, comm);
  for (int ptype = 0; ptype < swift_type_count; ++ptype)
    N_total[ptype] = offset[ptype] + N[ptype];

  /* The last rank now has the correct N_total. Let's
   * broadcast from there */
  MPI_Bcast(&N_total, 6, MPI_LONG_LONG_INT, mpi_size - 1, comm);

  /* Now everybody konws its offset and the total number of
   * particles of each type */

#ifdef IO_SPEED_MEASUREMENT
  ticks tic = getticks();
#endif

  /* Rank 0 prepares the file */
  if (mpi_rank == 0)
    prepare_velociraptor_snippet_file(e, baseName, N_total, ifile_mbp, internal_units, snapshot_units);

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef IO_SPEED_MEASUREMENT
  if (engine_rank == 0)
    message("Preparing file on rank 0 took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
#endif

  /* HDF5 File name */
  char fileName[FILENAME_BUFFER_SIZE];
  if (e->snapshot_int_time_label_on)
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%06i.hdf5", baseName,
             (int)round(e->time));
  else
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%04i.hdf5", baseName,
             e->snapshot_output_count);

  /* Now write the top-level cell structure */
  hid_t h_file_cells = 0, h_grp_cells = 0;
  if (mpi_rank == 0) {

    /* Open the snapshot on rank 0 */
    h_file_cells = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
    if (h_file_cells < 0)
      error("Error while opening file '%s' on rank %d.", fileName, mpi_rank);

    /* Create the group we want in the file */
    h_grp_cells = H5Gcreate(h_file_cells, "/Cells", H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (h_grp_cells < 0) error("Error while creating cells group");
  }

  /* Write the location of the particles in the arrays */
  io_write_cell_offsets(h_grp_cells, e->s->cdim, e->s->cells_top,
                        e->s->nr_cells, e->s->width, mpi_rank, N_total, offset,
                        internal_units, snapshot_units);

  /* Close everything */
  if (mpi_rank == 0) {
    H5Gclose(h_grp_cells);
    H5Fclose(h_file_cells);
  }

  /* Prepare some file-access properties */
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

  /* Set some MPI-IO parameters */
  // MPI_Info_set(info, "IBM_largeblock_io", "true");
  MPI_Info_set(info, "romio_cb_write", "enable");
  MPI_Info_set(info, "romio_ds_write", "disable");

  /* Activate parallel i/o */
  hid_t h_err = H5Pset_fapl_mpio(plist_id, comm, info);
  if (h_err < 0) error("Error setting parallel i/o");

  /* Align on 4k pages. */
  h_err = H5Pset_alignment(plist_id, 1024, 4096);
  if (h_err < 0) error("Error setting Hdf5 alignment");

  /* Disable meta-data cache eviction */
  H5AC_cache_config_t mdc_config;
  mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
  h_err = H5Pget_mdc_config(plist_id, &mdc_config);
  if (h_err < 0) error("Error getting the MDC config");

  mdc_config.evictions_enabled = 0; /* false */
  mdc_config.incr_mode = H5C_incr__off;
  mdc_config.decr_mode = H5C_decr__off;
  mdc_config.flash_incr_mode = H5C_flash_incr__off;
  h_err = H5Pset_mdc_config(plist_id, &mdc_config);
  if (h_err < 0) error("Error setting the MDC config");

/* Use parallel meta-data writes */
#if H5_VERSION_GE(1, 10, 0)
  h_err = H5Pset_all_coll_metadata_ops(plist_id, 1);
  if (h_err < 0) error("Error setting collective meta-data on all ops");
    // h_err = H5Pset_coll_metadata_write(plist_id, 1);
    // if (h_err < 0) error("Error setting collective meta-data writes");
#endif

#ifdef IO_SPEED_MEASUREMENT
  MPI_Barrier(MPI_COMM_WORLD);
  if (engine_rank == 0)
    message("Setting parallel HDF5 access properties took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
#endif

  /* Open HDF5 file with the chosen parameters */
  hid_t h_file = H5Fopen(fileName, H5F_ACC_RDWR, plist_id);
  if (h_file < 0) error("Error while opening file '%s'.", fileName);

#ifdef IO_SPEED_MEASUREMENT
  MPI_Barrier(MPI_COMM_WORLD);
  if (engine_rank == 0)
    message("Opening HDF5 file  took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
#endif

  /* Tell the user if a conversion will be needed */
  if (e->verbose && mpi_rank == 0) {
    if (units_are_equal(snapshot_units, internal_units)) {

      message("Snapshot and internal units match. No conversion needed.");

    } else {

      message("Conversion needed from:");
      message("(Snapshot) Unit system: U_M =      %e g.",
              snapshot_units->UnitMass_in_cgs);
      message("(Snapshot) Unit system: U_L =      %e cm.",
              snapshot_units->UnitLength_in_cgs);
      message("(Snapshot) Unit system: U_t =      %e s.",
              snapshot_units->UnitTime_in_cgs);
      message("(Snapshot) Unit system: U_I =      %e A.",
              snapshot_units->UnitCurrent_in_cgs);
      message("(Snapshot) Unit system: U_T =      %e K.",
              snapshot_units->UnitTemperature_in_cgs);
      message("to:");
      message("(internal) Unit system: U_M = %e g.",
              internal_units->UnitMass_in_cgs);
      message("(internal) Unit system: U_L = %e cm.",
              internal_units->UnitLength_in_cgs);
      message("(internal) Unit system: U_t = %e s.",
              internal_units->UnitTime_in_cgs);
      message("(internal) Unit system: U_I = %e A.",
              internal_units->UnitCurrent_in_cgs);
      message("(internal) Unit system: U_T = %e K.",
              internal_units->UnitTemperature_in_cgs);
    }
  }

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (N_total[ptype] == 0) continue;

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    hid_t h_grp = H5Gopen(h_file, partTypeGroupName, H5P_DEFAULT);
    if (h_grp < 0)
      error("Error while opening particle group %s.", partTypeGroupName);

    int num_fields = 0;
    struct io_props list[100];
    size_t Nparticles = 0;

    struct part* parts_written = NULL;
    struct xpart* xparts_written = NULL;
    struct gpart* gparts_written = NULL;
    struct spart* sparts_written = NULL;
    struct bpart* bparts_written = NULL;
    struct velociraptor_gpart_data* gpart_group_data_written = NULL;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas: {
          /* Ok, we need to fish out the particles we want */
          Nparticles = Ngas_written;

          /* Allocate temporary arrays */
          if (swift_memalign("parts_written", (void**)&parts_written,
                             part_align,
                             Ngas_written * sizeof(struct part)) != 0)
            error("Error while allocating temporary memory for parts");
          if (swift_memalign("xparts_written", (void**)&xparts_written,
                             xpart_align,
                             Ngas_written * sizeof(struct xpart)) != 0)
            error("Error while allocating temporary memory for xparts");

          /* Collect the particles we want to write */
          offset = parts->gpart - e->s->gparts;
          io_collect_velociraptor_parts_to_write(ifile_mbp, group_info+offset, parts, xparts, parts_written,
                                    xparts_written, Ngas, Ngas_written);

          /* Select the fields to write */
          hydro_write_particles(parts_written, xparts_written, list,
                                &num_fields);
          num_fields +=
              chemistry_write_particles(parts_written, list + num_fields);
          if (with_cooling || with_temperature) {
            num_fields +=
                cooling_write_particles(parts_written, xparts_written,
                                        list + num_fields, e->cooling_func);
          }
          num_fields += velociraptor_write_parts(
              parts_written, xparts_written, list + num_fields);
          num_fields += tracers_write_particles(
              parts_written, xparts_written, list + num_fields, with_cosmology);
          num_fields += star_formation_write_particles(
              parts_written, xparts_written, list + num_fields);
      } break;

      case swift_type_dark_matter: {
          Nparticles = Ndm_written;
          /* Allocate temporary array */
          if (swift_memalign("gparts_written", (void**)&gparts_written,
                             gpart_align,
                             Ndm_written * sizeof(struct gpart)) != 0)
            error("Error while allocating temporary memory for gparts");

            if (swift_memalign(
                    "gpart_group_written", (void**)&gpart_group_data_written,
                    gpart_align,
                    Ndm_written * sizeof(struct velociraptor_gpart_data)) != 0)
              error(
                  "Error while allocating temporary memory for gparts STF "
                  "data");

          /* Collect the non-inhibited DM particles from gpart */
          offset=0;
          io_collect_velociraptor_gparts_to_write(ifile_mbp, group_info+offset, gparts,
                                     gparts_written, gpart_group_data_written,
                                     Ntot, Ndm_written, with_stf);

          /* Select the fields to write */
          darkmatter_write_particles(ifile_mbp, gparts_written, list, &num_fields);
            num_fields += velociraptor_write_gparts(gpart_group_data_written,
                                                    list + num_fields);
      } break;

      case swift_type_stars: {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Nstars_written;

          /* Allocate temporary arrays */
          if (swift_memalign("sparts_written", (void**)&sparts_written,
                             spart_align,
                             Nstars_written * sizeof(struct spart)) != 0)
            error("Error while allocating temporary memory for sparts");

          /* Collect the particles we want to write */
          offset = sparts->gpart - e->s->gparts;
          io_collect_velociraptor_sparts_to_write(ifile_mbp, group_info+offset, sparts, sparts_written, Nstars,
                                     Nstars_written);

          /* Select the fields to write */
          stars_write_particles(sparts_written, list, &num_fields);
          num_fields += chemistry_write_sparticles(sparts, list + num_fields);
          num_fields += tracers_write_sparticles(sparts, list + num_fields,
                                                 with_cosmology);
          num_fields += fof_write_sparts(sparts_written, list + num_fields);
            num_fields +=
                velociraptor_write_sparts(sparts_written, list + num_fields);
      } break;

      case swift_type_black_hole: {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Nblackholes_written;

          /* Allocate temporary arrays */
          if (swift_memalign("bparts_written", (void**)&bparts_written,
                             bpart_align,
                             Nblackholes_written * sizeof(struct bpart)) != 0)
            error("Error while allocating temporary memory for bparts");

          /* Collect the particles we want to write */
          offset=bparts->gpart-e->s->gparts;
          io_collect_velociraptor_bparts_to_write(ifile_mbp, group_info+offset, bparts, bparts_written, Nblackholes,
                                     Nblackholes_written);

          /* Select the fields to write */
          black_holes_write_particles(bparts_written, list, &num_fields);
          num_fields += chemistry_write_bparticles(bparts, list + num_fields);
          num_fields += fof_write_bparts(bparts_written, list + num_fields);
            num_fields +=
                velociraptor_write_bparts(bparts_written, list + num_fields);
      } break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Write everything that is not cancelled */
    for (int i = 0; i < num_fields; ++i) {

      /* Did the user cancel this field? */
      char field[PARSER_MAX_LINE_SIZE];
      sprintf(field, "SelectOutput:%.*s_%s", FIELD_BUFFER_SIZE, list[i].name,
              part_type_names[ptype]);
      int should_write = parser_get_opt_param_int(params, field, 1);

      if (should_write)
        writeArray(e, h_grp, fileName, partTypeGroupName, list[i], Nparticles,
                   N_total[ptype], mpi_rank, offset[ptype], internal_units,
                   snapshot_units);
    }

    /* Free temporary array */
    if (parts_written) swift_free("parts_written", parts_written);
    if (xparts_written) swift_free("xparts_written", xparts_written);
    if (gparts_written) swift_free("gparts_written", gparts_written);
    if (gpart_group_data_written)
      swift_free("gpart_group_written", gpart_group_data_written);
    if (sparts_written) swift_free("sparts_written", sparts_written);
    if (bparts_written) swift_free("bparts_written", bparts_written);

#ifdef IO_SPEED_MEASUREMENT
    MPI_Barrier(MPI_COMM_WORLD);
    tic = getticks();
#endif

    /* Close particle group */
    H5Gclose(h_grp);

#ifdef IO_SPEED_MEASUREMENT
    MPI_Barrier(MPI_COMM_WORLD);
    if (engine_rank == 0)
      message("Closing particle group took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());

    tic = getticks();
#endif
  }

#ifdef IO_SPEED_MEASUREMENT
  MPI_Barrier(MPI_COMM_WORLD);
  tic = getticks();
#endif

  /* message("Done writing particles..."); */

  /* Close property descriptor */
  H5Pclose(plist_id);

#ifdef IO_SPEED_MEASUREMENT
  MPI_Barrier(MPI_COMM_WORLD);
  if (engine_rank == 0)
    message("Closing property descriptor took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();
#endif

  /* Close file */
  H5Fclose(h_file);

#ifdef IO_SPEED_MEASUREMENT
  MPI_Barrier(MPI_COMM_WORLD);
  if (engine_rank == 0)
    message("Closing file took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#endif

  e->snapshot_output_count++;
}

#endif

#if defined(HAVE_HDF5) && defined(WITH_MPI) && !defined(HAVE_PARALLEL_HDF5)

/**
* @brief Writes an serial HDF5 snippet output file (GADGET-3 type) with
*its XMF descriptor based on VR results.
 *
 * @param e The engine containing all the system.
 * @param baseName The common part of the snapshot file name.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 * @param mpi_rank The MPI rank of this node.
 * @param mpi_size The number of MPI ranks.
 * @param comm The MPI communicator.
 * @param info The MPI information object
 *
 * Creates an HDF5 output file and writes the particles contained
 * in the engine. If such a file already exists, it is erased and replaced
 * by the new one.
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output_serial(struct engine* e, const char* baseName,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units, int mpi_rank,
                         int mpi_size, MPI_Comm comm, MPI_Info info) {

  hid_t h_file = 0, h_grp = 0;
  int numFiles = 1;
  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;
  const struct gpart* gparts = e->s->gparts;
  const struct spart* sparts = e->s->sparts;
  const struct bpart* bparts = e->s->bparts;
  struct swift_params* params = e->parameter_file;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
  const int ifile_mbp = e->stf_output_mostboundparticles;
  ptrdiff_t offset;
  FILE* xmfFile = 0;

  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Nblackholes = e->s->nr_bparts;
  // const size_t Nbaryons = Ngas + Nstars;
  // const size_t Ndm = Ntot > 0 ? Ntot - Nbaryons : 0;

  /* Number of particles that we will write */
  size_t Ntot_written, Ngas_written, Nstars_written, Nblackholes_written;
  if (ifile_mbp) {
      for (int i=0; i<e->s->nr_gparts;i++) {
        if (group_info[i].groupID_lastmbp==0) continue;
        if (gparts[i].type==swift_type_gas) Ngas_written++;
        else if (gparts[i].type==swift_type_stars) Nstars_written++;
        else if (gparts[i].type==swift_type_black_hole) Nblackholes_written++;
        Ntot_written++;

      }
  }
  else {
      for (int i=0; i<e->s->nr_gparts;i++) {
        if (group_info[i].groupID==0) continue;
        if (gparts[i].type==swift_type_gas) Ngas_written++;
        else if (gparts[i].type==swift_type_stars) Nstars_written++;
        else if (gparts[i].type==swift_type_black_hole) Nblackholes_written++;
        Ntot_written++;
      }
  }
  const size_t Nbaryons_written =
      Ngas_written + Nstars_written + Nblackholes_written;
  const size_t Ndm_written =
      Ntot_written > 0 ? Ntot_written - Nbaryons_written : 0;

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  if (e->snapshot_int_time_label_on)
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%06i.hdf5", baseName,
             (int)round(e->time));
  else
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%04i.hdf5", baseName,
             e->snapshot_output_count);

  /* Compute offset in the file and total number of particles */
  size_t N[swift_type_count] = {Ngas_written,   Ndm_written,        0, 0,
                                Nstars_written, Nblackholes_written};
  long long N_total[swift_type_count] = {0};
  long long offset[swift_type_count] = {0};
  MPI_Exscan(&N, &offset, swift_type_count, MPI_LONG_LONG_INT, MPI_SUM, comm);
  for (int ptype = 0; ptype < swift_type_count; ++ptype)
    N_total[ptype] = offset[ptype] + N[ptype];

  /* The last rank now has the correct N_total. Let's broadcast from there */
  MPI_Bcast(&N_total, 6, MPI_LONG_LONG_INT, mpi_size - 1, comm);

  /* Now everybody konws its offset and the total number of particles of each
   * type */

  /* Do common stuff first */
  if (mpi_rank == 0) {

    /* First time, we need to create the XMF file */
    if (e->snapshot_output_count == 0) xmf_create_file(baseName);

    /* Prepare the XMF file for the new entry */
    xmfFile = xmf_prepare_file(baseName);

    /* Write the part corresponding to this specific output */
    xmf_write_outputheader(xmfFile, fileName, e->time);

    /* Open file */
    /* message("Opening file '%s'.", fileName); */
    h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) error("Error while opening file '%s'.", fileName);

    /* Open header to write simulation properties */
    /* message("Writing file header..."); */
    h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating file header\n");

    /* Convert basic output information to snapshot units */
    const double factor_time =
        units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_TIME);
    const double factor_length = units_conversion_factor(
        internal_units, snapshot_units, UNIT_CONV_LENGTH);
    const double dblTime = e->time * factor_time;
    const double dim[3] = {e->s->dim[0] * factor_length,
                           e->s->dim[1] * factor_length,
                           e->s->dim[2] * factor_length};

    /* Print the relevant information and print status */
    io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
    io_write_attribute(h_grp, "Time", DOUBLE, &dblTime, 1);
    const int dimension = (int)hydro_dimension;
    io_write_attribute(h_grp, "Dimension", INT, &dimension, 1);
    io_write_attribute(h_grp, "Redshift", DOUBLE, &e->cosmology->z, 1);
    io_write_attribute(h_grp, "Scale-factor", DOUBLE, &e->cosmology->a, 1);
    io_write_attribute_s(h_grp, "Code", "SWIFT");
    if (ifile_mbp) io_write_attribute_s(h_grp, "SnapshotFormat", "VR_MostBoundParticle_Snippet");
    else io_write_attribute_s(h_grp, "SnapshotFormat", "VR_InGroup_Snippet");
    time_t tm = time(NULL);
    io_write_attribute_s(h_grp, "Snapshot date", ctime(&tm));
    io_write_attribute_s(h_grp, "RunName", e->run_name);

    /* GADGET-2 legacy values */
    /* Number of particles of each type */
    unsigned int numParticles[swift_type_count] = {0};
    unsigned int numParticlesHighWord[swift_type_count] = {0};
    for (int ptype = 0; ptype < swift_type_count; ++ptype) {
      numParticles[ptype] = (unsigned int)N_total[ptype];
      numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);
    }
    io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, N_total,
                       swift_type_count);
    io_write_attribute(h_grp, "NumPart_Total", UINT, numParticles,
                       swift_type_count);
    io_write_attribute(h_grp, "NumPart_Total_HighWord", UINT,
                       numParticlesHighWord, swift_type_count);
    double MassTable[6] = {0., 0., 0., 0., 0., 0.};
    io_write_attribute(h_grp, "MassTable", DOUBLE, MassTable, swift_type_count);
    unsigned int flagEntropy[swift_type_count] = {0};
    flagEntropy[0] = writeEntropyFlag();
    io_write_attribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                       swift_type_count);
    io_write_attribute(h_grp, "NumFilesPerSnapshot", INT, &numFiles, 1);

    /* Close header */
    H5Gclose(h_grp);

    /* Print the code version */
    io_write_code_description(h_file);

    /* Print the run's policy */
    io_write_engine_policy(h_file, e);

    /* Print the SPH parameters */
    if (e->policy & engine_policy_hydro) {
      h_grp = H5Gcreate(h_file, "/HydroScheme", H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating SPH group");
      hydro_props_print_snapshot(h_grp, e->hydro_properties);
      hydro_write_flavour(h_grp);
      H5Gclose(h_grp);
    }

    /* Print the subgrid parameters */
    h_grp = H5Gcreate(h_file, "/SubgridScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating subgrid group");
    entropy_floor_write_flavour(h_grp);
    cooling_write_flavour(h_grp, e->cooling_func);
    chemistry_write_flavour(h_grp);
    tracers_write_flavour(h_grp);
    H5Gclose(h_grp);

    /* Print the gravity parameters */
    if (e->policy & engine_policy_self_gravity) {
      h_grp = H5Gcreate(h_file, "/GravityScheme", H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating gravity group");
      gravity_props_print_snapshot(h_grp, e->gravity_properties);
      H5Gclose(h_grp);
    }

    /* Print the stellar parameters */
    if (e->policy & engine_policy_stars) {
      h_grp = H5Gcreate(h_file, "/StarsScheme", H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating stars group");
      stars_props_print_snapshot(h_grp, e->stars_properties);
      H5Gclose(h_grp);
    }

    /* Print the cosmological model */
    h_grp =
        H5Gcreate(h_file, "/Cosmology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating cosmology group");
    if (e->policy & engine_policy_cosmology)
      io_write_attribute_i(h_grp, "Cosmological run", 1);
    else
      io_write_attribute_i(h_grp, "Cosmological run", 0);
    cosmology_write_model(h_grp, e->cosmology);
    H5Gclose(h_grp);

    /* Print the runtime parameters */
    h_grp =
        H5Gcreate(h_file, "/Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating parameters group");
    parser_write_params_to_hdf5(e->parameter_file, h_grp, 1);
    H5Gclose(h_grp);

    /* Print the runtime unused parameters */
    h_grp = H5Gcreate(h_file, "/UnusedParameters", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating parameters group");
    parser_write_params_to_hdf5(e->parameter_file, h_grp, 0);
    H5Gclose(h_grp);

    /* Print the system of Units used in the spashot */
    io_write_unit_system(h_file, snapshot_units, "Units");

    /* Print the system of Units used internally */
    io_write_unit_system(h_file, internal_units, "InternalCodeUnits");

    /* Tell the user if a conversion will be needed */
    if (e->verbose) {
      if (units_are_equal(snapshot_units, internal_units)) {

        message("Snapshot and internal units match. No conversion needed.");

      } else {

        message("Conversion needed from:");
        message("(Snapshot) Unit system: U_M =      %e g.",
                snapshot_units->UnitMass_in_cgs);
        message("(Snapshot) Unit system: U_L =      %e cm.",
                snapshot_units->UnitLength_in_cgs);
        message("(Snapshot) Unit system: U_t =      %e s.",
                snapshot_units->UnitTime_in_cgs);
        message("(Snapshot) Unit system: U_I =      %e A.",
                snapshot_units->UnitCurrent_in_cgs);
        message("(Snapshot) Unit system: U_T =      %e K.",
                snapshot_units->UnitTemperature_in_cgs);
        message("to:");
        message("(internal) Unit system: U_M = %e g.",
                internal_units->UnitMass_in_cgs);
        message("(internal) Unit system: U_L = %e cm.",
                internal_units->UnitLength_in_cgs);
        message("(internal) Unit system: U_t = %e s.",
                internal_units->UnitTime_in_cgs);
        message("(internal) Unit system: U_I = %e A.",
                internal_units->UnitCurrent_in_cgs);
        message("(internal) Unit system: U_T = %e K.",
                internal_units->UnitTemperature_in_cgs);
      }
    }

    /* Loop over all particle types */
    for (int ptype = 0; ptype < swift_type_count; ptype++) {

      /* Don't do anything if no particle of this kind */
      if (N_total[ptype] == 0) continue;

      /* Open the particle group in the file */
      char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
      snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
               ptype);
      h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (h_grp < 0) error("Error while creating particle group.\n");

      /* Close particle group */
      H5Gclose(h_grp);
    }

    /* Close file */
    H5Fclose(h_file);
  }

  /* Now write the top-level cell structure */
  hid_t h_file_cells = 0, h_grp_cells = 0;
  if (mpi_rank == 0) {

    /* Open the snapshot on rank 0 */
    h_file_cells = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
    if (h_file_cells < 0)
      error("Error while opening file '%s' on rank %d.", fileName, mpi_rank);

    /* Create the group we want in the file */
    h_grp_cells = H5Gcreate(h_file_cells, "/Cells", H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (h_grp_cells < 0) error("Error while creating cells group");
  }

  /* Write the location of the particles in the arrays */
  io_write_cell_offsets(h_grp_cells, e->s->cdim, e->s->cells_top,
                        e->s->nr_cells, e->s->width, mpi_rank, N_total, offset,
                        internal_units, snapshot_units);

  /* Close everything */
  if (mpi_rank == 0) {
    H5Gclose(h_grp_cells);
    H5Fclose(h_file_cells);
  }

  /* Now loop over ranks and write the data */
  for (int rank = 0; rank < mpi_size; ++rank) {

    /* Is it this rank's turn to write ? */
    if (rank == mpi_rank) {

      h_file = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
      if (h_file < 0)
        error("Error while opening file '%s' on rank %d.", fileName, mpi_rank);

      /* Loop over all particle types */
      for (int ptype = 0; ptype < swift_type_count; ptype++) {

        /* Don't do anything if no particle of this kind */
        if (N_total[ptype] == 0) continue;

        /* Add the global information for that particle type to the XMF
         * meta-file */
        if (mpi_rank == 0)
          xmf_write_groupheader(xmfFile, fileName, N_total[ptype],
                                (enum part_type)ptype);

        /* Open the particle group in the file */
        char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
        snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
                 ptype);
        h_grp = H5Gopen(h_file, partTypeGroupName, H5P_DEFAULT);
        if (h_grp < 0)
          error("Error while opening particle group %s.", partTypeGroupName);

        int num_fields = 0;
        struct io_props list[100];
        size_t Nparticles = 0;

        struct part* parts_written = NULL;
        struct xpart* xparts_written = NULL;
        struct gpart* gparts_written = NULL;
        struct spart* sparts_written = NULL;
        struct bpart* bparts_written = NULL;
        struct velociraptor_gpart_data* gpart_group_data_written = NULL;

        /* Write particle fields from the particle structure */
        switch (ptype) {

          case swift_type_gas: {
              Nparticles = Ngas_written;
              /* Allocate temporary arrays */
              if (swift_memalign("parts_written", (void**)&parts_written,
                                 part_align,
                                 Ngas_written * sizeof(struct part)) != 0)
                error("Error while allocating temporary memory for parts");
              if (swift_memalign("xparts_written", (void**)&xparts_written,
                                 xpart_align,
                                 Ngas_written * sizeof(struct xpart)) != 0)
                error("Error while allocating temporary memory for xparts");

              /* Collect the particles we want to write */
              offset = parts->gpart - e->s->gparts;
              io_collect_velociraptor_parts_to_write(ifile_mbp, group_info+offset, parts, xparts, parts_written,
                                        xparts_written, Ngas, Ngas_written);

              /* Select the fields to write */
              hydro_write_particles(parts_written, xparts_written, list,
                                    &num_fields);
              num_fields +=
                  chemistry_write_particles(parts_written, list + num_fields);
              if (with_cooling || with_temperature) {
                num_fields +=
                    cooling_write_particles(parts_written, xparts_written,
                                            list + num_fields, e->cooling_func);
              }
                num_fields += velociraptor_write_parts(
                    parts_written, xparts_written, list + num_fields);
              num_fields +=
                  tracers_write_particles(parts_written, xparts_written,
                                          list + num_fields, with_cosmology);
              num_fields += star_formation_write_particles(
                  parts_written, xparts_written, list + num_fields);
          } break;

          case swift_type_dark_matter: {
              /* Ok, we need to fish out the particles we want */
              Nparticles = Ndm_written;

              /* Allocate temporary array */
              if (swift_memalign("gparts_written", (void**)&gparts_written,
                                 gpart_align,
                                 Ndm_written * sizeof(struct gpart)) != 0)
                error("Error while allocating temporary memory for gparts");

                if (swift_memalign(
                        "gpart_group_written",
                        (void**)&gpart_group_data_written, gpart_align,
                        Ndm_written * sizeof(struct velociraptor_gpart_data)) !=
                    0)
                  error(
                      "Error while allocating temporary memory for gparts STF "
                      "data");

              /* Collect the non-inhibited DM particles from gpart */
              offset=0;
              io_collect_velociraptor_gparts_to_write(ifile_mbp, group_info+offset,
                  gparts, gparts_written,
                  gpart_group_data_written, Ntot, Ndm_written);

              /* Select the fields to write */
              darkmatter_write_particles(gparts_written, list, &num_fields);
                num_fields += velociraptor_write_gparts(
                    gpart_group_data_written, list + num_fields);
          } break;

          case swift_type_stars: {
              /* Ok, we need to fish out the particles we want */
              Nparticles = Nstars_written;

              /* Allocate temporary arrays */
              if (swift_memalign("sparts_written", (void**)&sparts_written,
                                 spart_align,
                                 Nstars_written * sizeof(struct spart)) != 0)
                error("Error while allocating temporary memory for sparts");

              /* Collect the particles we want to write */
              offset = sparts->gpart - e->s->gparts;
              io_collect_velociraptor_sparts_to_write(ifile_mbp, group_info+offset, sparts, sparts_written, Nstars,
                                         Nstars_written);

              /* Select the fields to write */
              stars_write_particles(sparts_written, list, &num_fields);
              num_fields +=
                  chemistry_write_sparticles(sparts_written, list + num_fields);
              num_fields += tracers_write_sparticles(
                  sparts_written, list + num_fields, with_cosmology);
                num_fields += velociraptor_write_sparts(sparts_written,
                                                        list + num_fields);
          } break;

          case swift_type_black_hole: {
              /* Ok, we need to fish out the particles we want */
              Nparticles = Nblackholes_written;

              /* Allocate temporary arrays */
              if (swift_memalign(
                      "bparts_written", (void**)&bparts_written, bpart_align,
                      Nblackholes_written * sizeof(struct bpart)) != 0)
                error("Error while allocating temporary memory for bparts");

              /* Collect the particles we want to write */
              offset = bparts->gpart - e->s->gparts;
              io_collect_velociraptor_bparts_to_write(ifile_mbp, group_info+offset, bparts, bparts_written, Nblackholes,
                                         Nblackholes_written);

              /* Select the fields to write */
              black_holes_write_particles(bparts_written, list, &num_fields);
              num_fields +=
                  chemistry_write_bparticles(bparts, list + num_fields);
                num_fields += velociraptor_write_bparts(bparts_written,
                                                        list + num_fields);
          } break;

          default:
            error("Particle Type %d not yet supported. Aborting", ptype);
        }

        /* Write everything that is not cancelled */
        for (int i = 0; i < num_fields; ++i) {

          /* Did the user cancel this field? */
          char field[PARSER_MAX_LINE_SIZE];
          sprintf(field, "SelectOutput:%.*s_%s", FIELD_BUFFER_SIZE,
                  list[i].name, part_type_names[ptype]);
          int should_write = parser_get_opt_param_int(params, field, 1);

          if (should_write)
            writeArray(e, h_grp, fileName, xmfFile, partTypeGroupName, list[i],
                       Nparticles, N_total[ptype], mpi_rank, offset[ptype],
                       internal_units, snapshot_units);
        }

        /* Free temporary array */
        if (parts_written) swift_free("parts_written", parts_written);
        if (xparts_written) swift_free("xparts_written", xparts_written);
        if (gparts_written) swift_free("gparts_written", gparts_written);
        if (gpart_group_data_written)
          swift_free("gpart_group_written", gpart_group_data_written);
        if (sparts_written) swift_free("sparts_written", sparts_written);
        if (bparts_written) swift_free("bparts_written", sparts_written);

        /* Close particle group */
        H5Gclose(h_grp);

        /* Close this particle group in the XMF file as well */
        if (mpi_rank == 0)
          xmf_write_groupfooter(xmfFile, (enum part_type)ptype);
      }

      /* Close file */
      H5Fclose(h_file);
    }

    /* Wait for the read of the reading to complete */
    MPI_Barrier(comm);
  }

  /* Write footer of LXMF file descriptor */
  if (mpi_rank == 0)
    xmf_write_outputfooter(xmfFile, e->snapshot_output_count, e->time);

  /* message("Done writing particles..."); */
  e->snapshot_output_count++;
}

#endif


#endif /* SWIFT_VELOCIRAPTOR_IO_H */
