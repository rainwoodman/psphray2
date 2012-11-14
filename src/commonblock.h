typedef struct {
    double Time;
    double a;
    double BoxSize;

    int SnapNumMajor;
    int SnapNumMajorBegin;
    int SnapNumMajorEnd;
    int SnapNumMinor;
    int IDByteSize;
    int NReader;

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

#define barrier() MPI_Barrier(MPI_COMM_WORLD)
#define abort(x) MPI_Abort(MPI_COMM_WORLD, (x))

#include "fckey.h"
#include "par.h"
#include "paramfile.h"
#include "snapshot.h"
