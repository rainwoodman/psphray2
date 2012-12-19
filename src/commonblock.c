#include <glib.h>
#include "mpiu.h"
#include "commonblock.h"

CommonBlock CB = {0};

void common_block_sync() {
    MPI_Bcast(&CB, sizeof(CB), MPI_BYTE, 0, MPI_COMM_WORLD);
    mpiu_bcast_string(&CB.datadir);
    mpiu_bcast_string(&CB.snapbase);
    mpiu_bcast_string(&CB.inputbase);
}
