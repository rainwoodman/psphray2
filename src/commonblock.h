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
    int SubDomainsPerTask;
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
    struct {
        double CritPhysDensity;
        double CritOverDensity;
        double FactorSN;
        double FactorEVP;
        double TempSupernova;
        double TempClouds;
        double MaxSfrTimescale;
    } SFREFF;
} CommonBlock;

extern CommonBlock CB;
void common_block_sync();

