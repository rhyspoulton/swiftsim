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

void io_collect_velociraptor_parts_to_write(const int imbp,
                                           const struct velociraptor_gpart_data* restrict groupinfo,
                                           const struct part* restrict parts,
                                           const struct xpart* restrict xparts,
                                           struct part* restrict parts_written,
                                           struct xpart* restrict xparts_written,
                                           const size_t Nparts,
                                           const size_t Nparts_written);
void io_collect_velociraptor_sparts_to_write(const int imbp,
                                            const struct velociraptor_gpart_data* restrict groupinfo,
                                            const struct spart* restrict sparts,
                                            struct spart* restrict sparts_written,
                                            const size_t Nsparts,
                                            const size_t Nsparts_written);
void io_collect_velociraptor_bparts_to_write(const int imbp,
                                            const struct velociraptor_gpart_data* restrict groupinfo,
                                            const struct bpart* restrict bparts,
                                            struct bpart* restrict bparts_written,
                                            const size_t Nbparts,
                                            const size_t Nbparts_written);
void io_collect_velociraptor_gparts_to_write(const int imbp,
                                            const struct velociraptor_gpart_data* restrict groupinfo,
                                            const struct gpart* restrict gparts,
                                            struct gpart* restrict gparts_written,
                                            const size_t Ngparts, const size_t Ngparts_written);

#if defined(HAVE_HDF5) && defined(WITH_MPI) && defined(HAVE_PARALLEL_HDF5)
//need for parallel hdf write
void prepare_velociraptor_snippet_file(struct engine* e, const char* baseName,
                  long long N_total[6],
                  const int ifile_mbp,
                  const struct unit_system* internal_units,
                  const struct unit_system* snapshot_units);

void write_velociraptor_snippet_output_parallel(struct engine* e, const char* baseName,
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units,
                           int mpi_rank, int mpi_size, MPI_Comm comm,
                           MPI_Info info);
#endif

#if defined(HAVE_HDF5) && defined(WITH_MPI) && !defined(HAVE_PARALLEL_HDF5)
void write_velociraptor_snippet_output_serial(struct engine* e, const char* baseName,
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units,
                           int mpi_rank, int mpi_size, MPI_Comm comm,
                           MPI_Info info);

                           void read_ic_serial(char* fileName, const struct unit_system* internal_units,
                                               double dim[3], struct part** parts, struct gpart** gparts,
                                               struct spart** sparts, struct bpart** bparts, size_t* Ngas,
                                               size_t* Ngparts, size_t* Nstars, size_t* Nblackholes,
                                               int* flag_entropy, int with_hydro, int with_gravity,
                                               int with_stars, int with_black_holes, int cleanup_h,
                                               int cleanup_sqrt_a, double h, double a, int mpi_rank,
                                               int mpi_size, MPI_Comm comm, MPI_Info info, int n_threads,
                                               int dry_run);

#endif

#endif /* SWIFT_VELOCIRAPTOR_IO_H */
