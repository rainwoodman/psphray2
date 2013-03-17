#include <glib.h>
#include <string.h>

#include "par.h"
#include "mpiu.h"

/*
 * collective call to do the domain decomposition.
 * */
static void domain_stage_one_distribute(PackedPar * pack);

void domain(PackedPar * pack) {
    domain_stage_one_distribute(pack);
}


static void count(PackedPar * pack, ipos_t (*bdry)[3], ptrdiff_t N[]) {
    int i;
    for(i = 0; i < NTask; i++) {
        N[i] = pstore_pack_searchsorted_left(pack, bdry[i]);
    }
}
static void bisect(ipos_t left[][3], 
        ipos_t middle[][3], 
        ipos_t right[][3]) {
    int i;
    for(i = 0; i < NTask; i++) {
        ipos_bisect(left[i], right[i], middle[i]);
    }
}
static int adjust(ptrdiff_t Ntot, 
        ptrdiff_t local_count[],
        ipos_t left[][3], ipos_t middle[][3], ipos_t right[][3]) {
    int i;
    ptrdiff_t global_count[NTask];
    ptrdiff_t tol = Ntot / NTask * 0.1;
    int finished = 0;
    
    MPI_Allreduce(local_count, global_count, 
        NTask, MPI_PTRDIFF, MPI_SUM, 
        MPI_COMM_WORLD);
    for(i = 0; i < NTask; i++) {
        if(ipos_compare(left[i], right[i]) == 0) {
            finished ++;
            continue;
        }
        ptrdiff_t expectation = (i + 1)* Ntot / NTask;
        if(global_count[i] < expectation) {
            ipos_immediate_next(middle[i], left[i]);
        } else {
            right[i][0] = middle[i][0];
            right[i][1] = middle[i][1];
            right[i][2] = middle[i][2];
        }
    }
    return finished != NTask;
}
static void domain_stage_one_distribute(PackedPar * pack) {
    ipos_t left[NTask][3];
    ipos_t right[NTask][3];
    ipos_t middle[NTask][3];
    ptrdiff_t local_count[NTask];
    ptrdiff_t global_count[NTask];
    ptrdiff_t Ntot;
    int i;
    for(i = 0; i < NTask; i++) {
        left[i][0] = left[i][1] = left[i][2] = 0;
        right[i][0] = right[i][1] = right[i][2] = IPOS_LIMIT - 1;
    }
    left[NTask - 1][0] = left[NTask - 1][1] = 
        left[NTask - 1][2] = IPOS_LIMIT - 1;

    MPI_Allreduce(&pack->size, &Ntot, 
            1, MPI_PTRDIFF, MPI_SUM, 
            MPI_COMM_WORLD);

    pstore_pack_sort(pack); 

    do {
        bisect(left, middle, right);
        count(pack, middle, local_count);
    } while(adjust(Ntot, local_count, left, middle, right));

    MPI_Barrier(MPI_COMM_WORLD);

    ROOTONLY {
        for(i = 0; i < NTask; i++) {
            char * is = ipos_str(left[i]);
            g_debug("%d, %s", i, is);
        }
    }
}
