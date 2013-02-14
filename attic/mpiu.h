#ifndef __MPIU_H__
#define __MPIU_H__
#include <mpi.h>

extern int ThisTask;
extern int NTask;
extern int PrevTask;
extern int NextTask;
#define ROOTONLY if(ThisTask == 0)
#define TAKETURNS \
    for(int __i__ = 0; MPI_Barrier(MPI_COMM_WORLD), __i__ < NTask; \
        MPI_Barrier(MPI_COMM_WORLD), __i__++) if(ThisTask == __i__)
void mpiu_module_init();
void mpiu_bcast_string(char ** string);
#endif
