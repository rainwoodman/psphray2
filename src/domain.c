#include <glib.h>
#include <mpi.h>

#include "commonblock.h"

static void cumbincount(fckey_t * POV, intptr_t * N) {
    for(int i = 0; i < NTask - 1; i++) {
        N[i] = par_search_by_fckey(&POV[i], PAR_BUFFER_IN);
    }
    N[NTask - 1] = NPARin;
    MPI_Allreduce(MPI_IN_PLACE, N, NTask, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
}
static void bincount_local(fckey_t * POV, intptr_t * N) {
    for(int i = 0; i < NTask - 1; i++) {
        N[i] = par_search_by_fckey(&POV[i], PAR_BUFFER_IN);
    }
    N[NTask - 1] = NPARin;
    for(int i = NTask - 1; i > 0; i--) {
        N[i] = N[i] - N[i - 1];
    }
}
static void find_pov(fckey_t * POV) {
    int i;
    fckey_t * POVleft =  g_new0(fckey_t, NTask);
    fckey_t * POVright = g_new0(fckey_t, NTask);
    fckey_t * POVmid = g_new0(fckey_t, NTask);

    for(i = 0; i < NTask; i++) {
        fckey_set_max(&POVright[i]);
        fckey_set_max(&POVmid[i]);
    }
    intptr_t * cumNtot_desired = g_new0(intptr_t, NTask);
    intptr_t * cumNtot = g_new0(intptr_t, NTask);

    intptr_t Ntot;
    MPI_Allreduce(&NPARin, &Ntot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    for(i = 0; i < NTask; i++) {
        cumNtot_desired[i] = (i + 1)* Ntot / NTask;
    }

    intptr_t tol = Ntot / NTask * CB.MemImbalanceTol;
    int k = 0;
    while(-1) {
        for(i = 0; i < NTask -1; i++) {
            fckey_center(&POVmid[i], &POVleft[i], &POVright[i]);
        }
        cumbincount(POVmid, cumNtot);
        for(i = 0; i < NTask; i++) {
            if(cumNtot[i] < cumNtot_desired[i] - tol) {
                POVleft[i] = POVmid[i];
                continue;
            }
            if(cumNtot[i] > cumNtot_desired[i] + tol) {
                POVright[i] = POVmid[i];
                continue;
            }
            POVleft[i] = POVmid[i];
            POVright[i] = POVmid[i];
        }
        int _break = TRUE;
        for(i = 0; i < NTask; i++) {
            if(fckey_cmp(&POVleft[i], &POVright[i])) {
                _break = FALSE;
            }
        }
        if(_break) break;
        k ++;
        if(k > 100) {
            ROOTONLY g_warning("too many iterations in domain decomposition");
            break;
        }
    }

    for(i = 0; i < NTask ; i++) {
        g_debug("task %d"
            #if 0
               " l " FCKEY_FMT " m " FCKEY_FMT " r " FCKEY_FMT
            #endif
               " dsr  %ld " 
               " real %ld " "\n",
               i, 
            #if 0
               FCKEY_PRINT(POVleft[i]), FCKEY_PRINT(POVmid[i]),
               FCKEY_PRINT(POVright[i]),
            #endif
               cumNtot_desired[i], 
               cumNtot[i]);
    }

    for(i = 0; i < NTask ; i++) {
        POV[i] = POVmid[i];
    }

    g_free(cumNtot);
    g_free(cumNtot_desired);
    g_free(POVleft);
    g_free(POVright);
    g_free(POVmid);
}


void domain_decompose() {
    /* after domain_decompose, PAR is available */
    fckey_t * POV = g_new0(fckey_t, NTask);
    intptr_t * segN = g_new0(intptr_t, NTask);
    intptr_t * N = g_new0(intptr_t, NTask);

    find_pov(POV);

    bincount_local(POV, N);

    MPI_Allreduce(N, segN, NTask, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    intptr_t Ntot;

    for(int i = 0; i < NTask; i++) {
        Ntot += segN[i];
    }

    intptr_t total = Ntot * CB.MemAllocFactor / NTask;
    intptr_t before = (total - segN[ThisTask]) / 2;

    /* so we have a few slots before and after, 
     * to adjusted domain boundaries against the tree */
    par_allocate(total - before, before);
    NPAR = segN[ThisTask];


    ROOTONLY {
        intptr_t segNmax=segN[0], segNmin=segN[0];
        for(int i = 0; i < NTask; i++) {
            segNmax = MAX(segNmax, segN[i]);
            segNmin = MIN(segNmax, segN[i]);
        }
        g_message("max load %ld, min load %ld", segNmax, segNmin);
    }

    /* now lets decide the communication layout */
    
    int *sendcounts = g_new0(int, NTask);
    int *sdispls = g_new0(int, NTask);
    int *recvcounts = g_new0(int, NTask);
    int *rdispls = g_new0(int, NTask);

    for(int i = 0; i < NTask; i++) {
        sendcounts[i] = N[i];
    }

    MPI_Alltoall(sendcounts, 1, MPI_INT,
                 recvcounts, 1, MPI_INT,
                      MPI_COMM_WORLD);

    for(int i = 0; i < NTask; i++) {
        sdispls[i] = (i>0)?(sdispls[i-1] + sendcounts[i-1]):0;
        rdispls[i] = (i>0)?(rdispls[i-1] + recvcounts[i-1]):0;
    }


    MPI_Datatype partype = 0;
    MPI_Type_contiguous(sizeof(par_t), MPI_BYTE, &partype);
    MPI_Type_commit(&partype);

    MPI_Alltoallv(&PARin(0), sendcounts, sdispls, partype,
                &PAR(0), recvcounts, rdispls, partype,
                MPI_COMM_WORLD);

    MPI_Type_free(&partype);
    g_free(rdispls);
    g_free(recvcounts);
    g_free(sdispls);
    g_free(sendcounts);
    g_free(segN);
    g_free(POV);
}

void domain_adjust() {
    Node * first = tree_locate_fckey(TREEROOT, &PAR(0).fckey);
    Node * last = tree_locate_fckey(TREEROOT, &PAR(-1).fckey);
    Node before;
    Node behind;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Sendrecv(first, sizeof(Node), MPI_BYTE, (ThisTask - 1 + NTask) % NTask, 1,
        &behind, sizeof(Node), MPI_BYTE, (ThisTask + 1) % NTask, 1,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Sendrecv(last, sizeof(Node), MPI_BYTE, (ThisTask + 1) % NTask, 2,
        &before, sizeof(Node), MPI_BYTE, (ThisTask - 1 + NTask) % NTask, 2,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    intptr_t in_before = 0;
    intptr_t in_behind = 0;
    intptr_t from_before = 0;
    intptr_t from_behind = 0;

    for(intptr_t i = 0; i < NPAR; i++) {
        if(!tree_node_contains_fckey(&before, &PAR(i).fckey))
            break;
        in_before ++;
    }
    for(intptr_t i = NPAR - 1; i >= 0; i--) {
        if(!tree_node_contains_fckey(&behind, &PAR(i).fckey))
            break;
        in_behind ++;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Sendrecv(&in_before, 1, MPI_LONG, (ThisTask - 1 + NTask) % NTask, 1,
        &from_behind, 1, MPI_LONG, (ThisTask + 1) % NTask, 1,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Sendrecv(&in_behind, 1, MPI_LONG, (ThisTask + 1) % NTask, 2,
        &from_before, 1, MPI_LONG, (ThisTask - 1 + NTask) % NTask, 2,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    g_print("in_before %ld from_before %ld " NODE_FMT "\n", 
            in_before, from_before, NODE_PRINT(before));
    g_print("in_behind %ld from_behind %ld " NODE_FMT "\n", 
            in_behind, from_behind, NODE_PRINT(behind));
}
