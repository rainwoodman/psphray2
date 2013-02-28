#include <glib.h>
#include "mpiu.h"
#include "commonblock.h"

CommonBlock CB = {0};

void common_block_sync() {
    MPI_Bcast(&CB, sizeof(CB), MPI_BYTE, 0, MPI_COMM_WORLD);
    mpiu_bcast_string(&CB.datadir);
    ROOTONLY {
    } else {
        CB.IC.R = g_new0(region_t, CB.IC.NRegions);
        CB.IC.Levels = g_new0(level_t, CB.IC.NLevels);
    }
    MPI_Bcast(CB.IC.R, sizeof(region_t) * CB.IC.NRegions, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(CB.IC.Levels, sizeof(level_t) * CB.IC.NLevels, MPI_BYTE, 0, MPI_COMM_WORLD);
}
