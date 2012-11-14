#include <glib.h>
#include <mpi.h>
#include <string.h>
#include "commonblock.h"

CommonBlock CB = {0};
int ThisTask;
int NTask;

void bcast_string(char ** string) {
    int len;
    ROOTONLY {
        len = strlen(string[0]);
    }

    MPI_Bcast(&len, sizeof(len), MPI_BYTE, 0, MPI_COMM_WORLD);
    ROOTONLY { } else {
        string[0] = g_new(char, len + 1);
    }
    MPI_Bcast(string[0], len, MPI_BYTE, 0, MPI_COMM_WORLD);
}

void common_block_bootstrap() {
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
}

void common_block_sync() {
    MPI_Bcast(&CB, sizeof(CB), MPI_BYTE, 0, MPI_COMM_WORLD);
    bcast_string(&CB.datadir);
    bcast_string(&CB.snapbase);
    bcast_string(&CB.inputbase);
}
