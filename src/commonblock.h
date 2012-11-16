typedef struct {
    double Time;
    double a;

    double BoxSize; /* filled by snapshot.c, not by paramfile.c, yet */

    int SnapNumMajor;
    int SnapNumMajorBegin;
    int SnapNumMajorEnd;
    int SnapNumMinor;
    int IDByteSize;
    int NReader;

    double MemImbalanceTol;
    double MemAllocFactor;

    int NodeSplitThresh;

    /* do not forget to add the bcast call to common_block_sync() */
    char * datadir;
    char * inputbase;
    char * snapbase;

    struct {
        int VERBOSE;
    } F;
    struct {
        double MYEAR_h;
        double KPC_h;

        double CM_h;
        double GRAM_h;
        double SECOND_h;

        double CM;
        double GRAM;
        double SECOND;
        double SOLARMASS;
        double PROTONMASS;
    } U;
    struct {
        double h;
        double H;
        double G;
        double C;
        double OmegaB;
        double OmegaM;
        double OmegaL;
    } C;
} CommonBlock;

extern CommonBlock CB;
extern int ThisTask;
extern int NTask;
void common_block_sync();
void common_block_bootstrap();

#define ROOTONLY if(ThisTask == 0)
#define TAKETURNS \
    for(int __i__ = 0; __i__ < NTask; \
        MPI_Barrier(MPI_COMM_WORLD), __i__++) if(ThisTask == __i__)


#include "fckey.h"
#include "par.h"
#include "paramfile.h"
#include "snapshot.h"
#include "tree.h"

