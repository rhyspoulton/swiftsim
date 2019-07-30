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
#ifndef SWIFT_VELOCIRAPTOR_STRUCT_H
#define SWIFT_VELOCIRAPTOR_STRUCT_H

/* Config parameters. */
#include "../config.h"

/**
 * @brief Data returned by VELOCIraptor for each #gpart.
 */
struct velociraptor_gpart_data {
  /* store the group id for the gravity particle */
  long long groupID;
  /* store whether the particle is the most bound particle of a group */
  bool imbp;
  /* store the groupID of the last group to which the particle was the group's most bound particle */
  long long groupID_lastmbp;
};

#endif /* SWIFT_VELOCIRAPTOR_STRUCT_H */
