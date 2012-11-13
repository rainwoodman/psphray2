typedef struct {
    int ThisTask;
    int NTask;
    double Time;
    double a;

    char * datadir;

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

#define ROOTONLY if(CB.ThisTask == 0)

#define barrier() MPI_Barrier(MPI_COMM_WORLD)
#define abort(x) MPI_Abort(MPI_COMM_WORLD, (x))

