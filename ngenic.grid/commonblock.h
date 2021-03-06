typedef struct {
    double center[3];
    double size[3];
    int ibottom[3];
    int itop[3];
    int isize[3];
    ptrdiff_t stride[3];
} region_t;

typedef struct {
    int Nmesh;
    double Scale;
    int DownSample;
} level_t;

typedef struct pow_table {
  double logk, logD;
} pow_table;


typedef struct {
    double a;

    double BoxSize; /* filled by snapshot.c, not by paramfile.c, yet */

    /* do not forget to add the bcast call to common_block_sync() */
    char * datadir;

    struct {
        int VERBOSE;
        int INDEX;
        int POWERONLY;
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
        int SphereMode;
        int WhichSpectrum;
        int NmeshPrimary;
        char * PowerSpectrumFile;
        double PrimordialIndex;
        double ShapeGamma;
    } IC;
} CommonBlock;

extern CommonBlock CB;
void common_block_sync();

extern int NR;
extern int NL;
int NPowerTable;
extern region_t * R;
extern level_t * L;
extern pow_table *PowerTable;
void levels_sort();
void levels_alloc(size_t n);
void regions_alloc(size_t n);
void pow_table_alloc(size_t n);
int levels_select(int Nmesh);
