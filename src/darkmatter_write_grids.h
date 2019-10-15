#ifndef SWIFT_DARKMATTER_WRITE_GRIDS_H
#define SWIFT_DARKMATTER_WRITE_GRIDS_H

/* Config parameters. */
#include "../config.h"
#include "engine.h"
#include "units.h"

void darkmatter_write_grids(struct engine* e, const size_t Npart,
                            const hid_t h_file,
                            const struct unit_system* internal_units,
                            const struct unit_system* snapshot_units);

#endif /* SWIFT_DARKMATTER_WRITE_GRIDS_H */
