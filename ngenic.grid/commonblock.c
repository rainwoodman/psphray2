#include <glib.h>
#include "mpiu.h"
#include "commonblock.h"
#include <stdlib.h>

CommonBlock CB = {0};
int NR;
int NL;
int NPowerTable;
region_t * R;
level_t * L;
pow_table *PowerTable;

void common_block_sync() {
    MPI_Bcast(&CB, sizeof(CB), MPI_BYTE, 0, MPI_COMM_WORLD);
    mpiu_bcast_string(&CB.datadir);
    MPI_Bcast(&NR, sizeof(NR), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NL, sizeof(NL), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NPowerTable, sizeof(NPowerTable), MPI_BYTE, 0, MPI_COMM_WORLD);
    ROOTONLY {
    } else {
        R = g_new0(region_t, NR);
        L = g_new0(level_t, NL);
        PowerTable = g_new0(pow_table, NPowerTable);
    }
    MPI_Bcast(R, sizeof(region_t) * NR, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(L, sizeof(level_t) * NL, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(PowerTable, sizeof(pow_table) * NPowerTable, MPI_BYTE, 0, MPI_COMM_WORLD);
}

void regions_alloc(size_t n) {
    NR = n;
    R = g_new0(region_t, n);
}
void levels_alloc(size_t n) {
    NL = n;
    L = g_new0(level_t, n);
}
void pow_table_alloc(size_t n) {
    NPowerTable = n;
    PowerTable = g_new0(pow_table, n);
}

static int cmp_level_t_nmesh(level_t * l1, level_t *  l2) {
    return (l1->Nmesh > l2->Nmesh) - (l2->Nmesh > l1->Nmesh);
}
void levels_sort() {
    qsort(L, NL, sizeof(level_t), (int(*)(const void *, const void *) )cmp_level_t_nmesh);
}

int levels_select(int Nmesh) {
    int i;
    for(i = 0; i < NL; i++) {
        if (L[i].Nmesh == Nmesh) {
            return i;
        }
    }
    g_assert_not_reached();
}
