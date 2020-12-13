/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_CELL_STARS_H
#define SWIFT_CELL_STARS_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "lock.h"
#include "star_formation_logger_struct.h"
#include "timeline.h"

/**
 * @brief Stars-related cell variables.
 */
struct cell_stars {

  /*! Pointer to the #spart data. */
  struct spart *parts;

  /*! Pointer to the #spart data at rebuild time. */
  struct spart *parts_rebuild;

  /*! The star ghost task itself */
  struct task *ghost;

  /*! Linked list of the tasks computing this cell's star density. */
  struct link *density;

  /*! Linked list of the tasks computing this cell's star feedback. */
  struct link *feedback;

  /*! The task computing this cell's sorts before the density. */
  struct task *sorts;

  /*! The drift task for sparts */
  struct task *drift;

  /*! Implicit tasks marking the entry of the stellar physics block of tasks
   */
  struct task *stars_in;

  /*! Implicit tasks marking the exit of the stellar physics block of tasks */
  struct task *stars_out;

  /*! Last (integer) time the cell's spart were drifted forward in time. */
  integertime_t ti_old_part;

  /*! Spin lock for various uses (#spart case). */
  swift_lock_type lock;

  /*! Spin lock for star formation use. */
  swift_lock_type star_formation_lock;

  /*! Nr of #spart in this cell. */
  int count;

  /*! Nr of #spart this cell can hold after addition of new #spart. */
  int count_total;

  /*! Max smoothing length in this cell. */
  float h_max;

  /*! Values of h_max before the drifts, used for sub-cell tasks. */
  float h_max_old;

  /*! Maximum part movement in this cell since last construction. */
  float dx_max_part;

  /*! Values of dx_max before the drifts, used for sub-cell tasks. */
  float dx_max_part_old;

  /*! Maximum particle movement in this cell since the last sort. */
  float dx_max_sort;

  /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
  float dx_max_sort_old;

  /*! Pointer for the sorted indices. */
  struct sort_entry *sort;

  /*! Bit mask of sort directions that will be needed in the next timestep. */
  uint16_t requires_sorts;

  /*! Bit-mask indicating the sorted directions */
  uint16_t sorted;

  /*! Bit-mask indicating the sorted directions */
  uint16_t sort_allocated;

  /*! Bit mask of sorts that need to be computed for this cell. */
  uint16_t do_sort;

  /*! Maximum end of (integer) time step in this cell for star tasks. */
  integertime_t ti_end_min;

  /*! Maximum end of (integer) time step in this cell for star tasks. */
  integertime_t ti_end_max;

  /*! Maximum beginning of (integer) time step in this cell for star tasks.
   */
  integertime_t ti_beg_max;

  /*! Number of #spart updated in this cell. */
  int updated;

  /*! Is the #spart data of this cell being used in a sub-cell? */
  int hold;

  /*! Star formation history struct */
  struct star_formation_history sfh;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Last (integer) time the cell's sort arrays were updated. */
  integertime_t ti_sort;
#endif
};

#endif /* SWIFT_CELL_STARS_H */
