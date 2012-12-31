typedef struct {
    double center[3];
    double size[3];
} region_t;

typedef struct {
    double a;

    double BoxSize; /* filled by snapshot.c, not by paramfile.c, yet */

    /* do not forget to add the bcast call to common_block_sync() */
    char * datadir;

    struct {
        int VERBOSE;
        int INDEX;
        int FULL;
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
        double Sigma8;
    } C;
    struct {
        int Seed;
        int Nmesh;
        int WhichSpectrum;
        double PrimordialIndex;
        double ShapeGamma;
        double Scale;
        int NRegions;
        region_t * R;
    } IC;
} CommonBlock;

extern CommonBlock CB;
void common_block_sync();

