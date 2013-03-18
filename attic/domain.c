#include <glib.h>
#include <string.h>

#include "par.h"
#include "mpiu.h"

/*
 * collective call to do the domain decomposition.
 * */
static void domain_stage_one_distribute(PackedPar ** pack);

void domain(PackedPar * pack) {
    domain_stage_one_distribute(&pack);

    PStore * pstore = pstore_new(128);
    pstore_insert_pack(pstore, pack);
    pstore_pack_free(pack);
    MPI_Barrier(MPI_COMM_WORLD);
    g_message("%03d build tree", ThisTask);
}

/**
 * stage one will evenly distribute the particles
 * to all processes.
 * the particles are first z-ordered, then
 * the curve is divided into NTask pieces, and each
 * process collects a piece.
 * this is done by assuming not too many particles have
 * identical ipos.  (which would be a serious problem anyways)
 */
static PackedPar * exchange(PackedPar * pack, ptrdiff_t local_count[]);
static void count(PackedPar * pack, ipos_t (*bdry)[3], ptrdiff_t N[]);
static void bisect(ipos_t left[][3], 
        ipos_t middle[][3], 
        ipos_t right[][3]);
static int adjust(ptrdiff_t Ntot, 
        ptrdiff_t local_count[],
        ipos_t left[][3], ipos_t middle[][3], ipos_t right[][3]);

static void domain_stage_one_distribute(PackedPar ** ppack) {
    PackedPar * pack = *ppack;

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

    *ppack = exchange(pack, local_count);
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

#define MPI_TAG_NBYTES 0
#define MPI_TAG_PARBYTES 1
#define MPI_TAG_SIZE 2
#define MPI_TAG_PACK 3

PackedPar * exchange(PackedPar * pack, ptrdiff_t local_count[]) {
    int i;
    MPI_Request requests[4][NTask];
    PackedPar * sendbuf[NTask];
    ptrdiff_t recvcount[NTask];
    ptrdiff_t recvparbytes[NTask];
    ptrdiff_t recvnbytes[NTask];
    ptrdiff_t sum_recvcount = 0;
    ptrdiff_t sum_recvparbytes = 0;
    ptrdiff_t zero = 0; 

    pstore_pack_leak_check_start();
    for(i = 0; i < NTask; i++) {
        ptrdiff_t start = i?local_count[i - 1]: 0;
        ptrdiff_t end = local_count[i];
        ptrdiff_t size = end - start;
        if(end == start) {
            /* as send is non-blocking, the buffer has
             * to be avail before wait returns,
             * thus we feed a zero variable decleared before
             * the loop */
            MPI_Issend(&zero, 1, MPI_PTRDIFF, 
                        i, MPI_TAG_SIZE, MPI_COMM_WORLD,
                        &requests[0][i]);
            sendbuf[i] = NULL;
            requests[1][i] = MPI_REQUEST_NULL;
            requests[2][i] = MPI_REQUEST_NULL;
            requests[3][i] = MPI_REQUEST_NULL;
            continue;
        }
        sendbuf[i] = pstore_pack_create_from_selection(
                pack, start, end);
        MPI_Issend(&sendbuf[i]->size, 1, MPI_PTRDIFF, 
                    i, MPI_TAG_SIZE, MPI_COMM_WORLD,
                    &requests[0][i]);
        MPI_Issend(&sendbuf[i]->parbytes, 1, MPI_PTRDIFF, 
                    i, MPI_TAG_PARBYTES, MPI_COMM_WORLD, 
                    &requests[1][i]);
        MPI_Issend(&sendbuf[i]->nbytes, 1, MPI_PTRDIFF, 
                    i, MPI_TAG_NBYTES, MPI_COMM_WORLD, 
                    &requests[2][i]);
        MPI_Issend(sendbuf[i], sendbuf[i]->nbytes, MPI_BYTE, 
                    i, MPI_TAG_PACK, MPI_COMM_WORLD, 
                    &requests[3][i]);
    }
    
    for(i = 0; i < NTask; i++) {
        MPI_Recv(&recvcount[i], 1, MPI_PTRDIFF, i, MPI_TAG_SIZE, 
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(recvcount[i] > 0) {
            MPI_Recv(&recvparbytes[i], 1, MPI_PTRDIFF, i, MPI_TAG_PARBYTES, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&recvnbytes[i], 1, MPI_PTRDIFF, i, MPI_TAG_NBYTES, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            recvnbytes[i] = 0;
            recvparbytes[i] = 0;
        }
        sum_recvcount += recvcount[i];
        sum_recvparbytes += recvparbytes[i];
    }
    
    pstore_pack_free(pack);

    pack = pstore_pack_create_direct(sum_recvcount, sum_recvparbytes);
    ptrdiff_t cursor = 0;
    for(i = 0; i < NTask; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(recvcount[i] == 0) continue;
        PackedPar * recvbuf = pstore_pack_create_direct(recvcount[i],
                        recvparbytes[i]);
        g_assert(recvnbytes[i] == recvbuf->nbytes);
        MPI_Recv(recvbuf, recvnbytes[i], MPI_BYTE, 
                i, MPI_TAG_PACK, 
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
        ptrdiff_t j;
        for(j = 0; j < recvbuf->size; j++) {
            pstore_pack_push(pack, &cursor, pstore_pack_get(recvbuf, j));
        }
        pstore_pack_free(recvbuf);
    }

    MPI_Waitall(NTask, requests[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(NTask, requests[1], MPI_STATUSES_IGNORE);
    MPI_Waitall(NTask, requests[2], MPI_STATUSES_IGNORE);
    int index;
    while(MPI_UNDEFINED != (
            MPI_Waitany(NTask, requests[3], 
                &index, MPI_STATUS_IGNORE), 
            index)) {
        pstore_pack_free(sendbuf[index]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    g_message("%03d all send finished", ThisTask);
    pstore_pack_leak_check_end();
    return pack; 
}
