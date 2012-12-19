#include <glib.h>
#include <mpi.h>
#include <string.h>
#include "mpiu.h"
int ThisTask;
int NTask;
int PrevTask, NextTask;

void mpiu_init() {
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    PrevTask = (ThisTask - 1 + NTask) % NTask,
    NextTask = (ThisTask + 1) % NTask;
}

void mpiu_bcast_string(char ** string) {
    int len;
    ROOTONLY {
        len = strlen(string[0]);
    }

    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    ROOTONLY { } else {
        string[0] = g_new(char, len + 1);
    }
    MPI_Bcast(string[0], len + 1, MPI_BYTE, 0, MPI_COMM_WORLD);
}

