#include <glib.h>
#include <mpi.h>

#include "commonblock.h"

static void cumbincount(fckey_t * POV, intptr_t * N) {
    for(int i = 0; i < NTask - 1; i++) {
        N[i] = par_search_by_fckey(PAR_BUFFER_IN, &POV[i]);
    }
    N[NTask - 1] = NPARin;
    MPI_Allreduce(MPI_IN_PLACE, N, NTask, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
}
static void bincount_local(fckey_t * POV, intptr_t * N) {
    for(int i = 0; i < NTask - 1; i++) {
        N[i] = par_search_by_fckey(PAR_BUFFER_IN, &POV[i]);
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
        if(k > 500) {
            ROOTONLY g_warning("too many iterations in domain decomposition");
            break;
        }
    }

    ROOTONLY {
        for(i = 0; i < NTask ; i++) {
            g_message("task %d"
                #if 1
                   " l " FCKEY_FMT " m " FCKEY_FMT " r " FCKEY_FMT
                #endif
                   " dsr  %ld " 
                   " real %ld " "\n",
                   i, 
                #if 1
                   FCKEY_PRINT(POVleft[i]), FCKEY_PRINT(POVmid[i]),
                   FCKEY_PRINT(POVright[i]),
                #endif
                   cumNtot_desired[i], 
                   cumNtot[i]);
        }
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
    par_sort_by_fckey(PAR_BUFFER_IN);
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
    par_reserve(PAR_BUFFER_MAIN, total - before, before);
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

    par_sort_by_fckey(PAR_BUFFER_MAIN);
    par_free(PAR_BUFFER_IN);
}

static void domain_mark_complete() {
    fckey_t * before = g_slice_new(fckey_t);
    fckey_t * behind = g_slice_new(fckey_t);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Sendrecv(&PAR(0).fckey, sizeof(fckey_t), MPI_BYTE, PrevTask, 101,
        behind, sizeof(fckey_t), MPI_BYTE, NextTask, 101,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Sendrecv(&PAR(-1).fckey, sizeof(fckey_t), MPI_BYTE, NextTask, 102,
        before, sizeof(fckey_t), MPI_BYTE, PrevTask, 102,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);

    TreeIter iter;
    tree_iter_init(&iter, TREEROOT);
    intptr_t complete_count = 0;
    intptr_t incomplete_count = 0;
    for(Node * node = tree_iter_next(&iter);
        node;
        node = tree_iter_next(&iter)) {
        if(ThisTask != 0 && 
            tree_node_contains_fckey(node, before)) {
            node->complete = FALSE;
            incomplete_count ++;
          #if 0
            g_print("incomplete" NODE_FMT ", key" FCKEY_FMT "\n",
                NODE_PRINT(node[0]), FCKEY_PRINT(before[0]));
          #endif
            continue;
        }
        if(ThisTask != NTask -1 &&
            tree_node_contains_fckey(node, behind)) {
            node->complete = FALSE;
            incomplete_count ++;
          #if 0
            g_print("incomplete" NODE_FMT ", key" FCKEY_FMT "\n",
                NODE_PRINT(node[0]), FCKEY_PRINT(behind[0]));
          #endif
            continue;
        }
        node->complete = TRUE;
        complete_count ++;
    }
    g_slice_free(fckey_t, behind);
    g_slice_free(fckey_t, before);

    TAKETURNS {
        g_print("%02d complete %ld incomplete %ld\n", ThisTask,
            complete_count, incomplete_count);
    }
}

void domain_build_tree() {
    /* 
     *
     * adjust the domain boundary to build a tree, 
     * so that the tree nodes
     * do not lie across domains.
     *
     * The algorithm:
     *
     * Consider two adjacent trees, and their nearby edge nodes A B.
     * Two cases:
     *  a) A B are disjoint:
     *      no need to do anything.
     *  b) A B are joint (identical or children):
     *      assume B is lower than A.
     *      on B's task find the image of A.
     *      compare len(A') and len(A). migrate
     *      the shorter one.
     *  */

    tree_build();

    Node * first = tree_locate_fckey(TREEROOT, &PAR(0).fckey);
    Node * last = tree_locate_fckey(TREEROOT, &PAR(-1).fckey);
    Node * before = g_slice_new(Node);
    Node * behind = g_slice_new(Node);

    g_message("root" NODE_FMT "", NODE_PRINT(TREEROOT[0]));
    g_message("first" NODE_FMT "", NODE_PRINT(first[0]));
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Sendrecv(first, sizeof(Node), MPI_BYTE, PrevTask, 1,
        behind, sizeof(Node), MPI_BYTE, NextTask, 1,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Sendrecv(last, sizeof(Node), MPI_BYTE, NextTask, 2,
        before, sizeof(Node), MPI_BYTE, PrevTask, 2,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if(tree_node_contains_node(behind, last)) {
        last = tree_node_find_image(last, behind);
    }
    if(tree_node_contains_node(before, first)) {
        first = tree_node_find_image(first, before);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Sendrecv(first, sizeof(Node), MPI_BYTE, PrevTask, 3,
        behind, sizeof(Node), MPI_BYTE, NextTask, 3,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Sendrecv(last, sizeof(Node), MPI_BYTE, NextTask, 4,
        before, sizeof(Node), MPI_BYTE, PrevTask, 4,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Datatype partype;
    MPI_Type_contiguous(sizeof(par_t), MPI_BYTE, &partype);
    MPI_Type_commit(&partype);
    MPI_Request head_req;
    MPI_Request tail_req;
    int head_sendcount = 0;
    int head_recvcount = 0;
    int tail_sendcount = 0;
    int tail_recvcount = 0;

    int do_tail = tree_node_contains_node(last, behind) || 
                tree_node_contains_node(behind, last);
    int do_head = tree_node_contains_node(first, before) || 
                tree_node_contains_node(before, first);

    if(do_tail) {
        if(last->npar < behind->npar) {
            /* send */
            tail_sendcount = last->npar;
            par_t * sendbuf = par_append(PAR_BUFFER_MAIN, 
                                   -tail_sendcount);
            MPI_Isend(sendbuf, tail_sendcount, partype,
                NextTask, 6,
                MPI_COMM_WORLD, &tail_req);
        } else {
            /* recv */ 
            tail_recvcount = behind->npar;
            par_t * recvbuf = par_append(PAR_BUFFER_MAIN, 
                              tail_recvcount);
            MPI_Irecv(recvbuf, tail_recvcount, partype, 
                NextTask, 7, 
                MPI_COMM_WORLD, &tail_req);
        }
    }
    if(do_head) {
        if(first->npar > before->npar) {
            /* recv */
            head_recvcount = before->npar;
            par_t * recvbuf = par_prepend(PAR_BUFFER_MAIN, 
                                head_recvcount);
            MPI_Irecv(recvbuf, head_recvcount, partype, 
                PrevTask, 6,
                MPI_COMM_WORLD, &head_req);
        } else {
            /* send */
            head_sendcount = first->npar;
            par_t * sendbuf = par_prepend(PAR_BUFFER_MAIN, 
                               -head_sendcount);
            MPI_Isend(sendbuf, head_sendcount, partype, 
                PrevTask, 7,
                MPI_COMM_WORLD, &head_req);
        }
    }
    intptr_t exchange = head_sendcount + tail_sendcount;
    intptr_t total_exchange = 0;
    
    MPI_Allreduce(&exchange, &total_exchange, 1, MPI_LONG, 
            MPI_SUM, MPI_COMM_WORLD);

    ROOTONLY {
        g_message("exchange of %ld particles ", total_exchange);
    }
    if(total_exchange > 0) {
        TAKETURNS {
            g_print("%02d -> %02d (%05d) "
                    "%02d <- %02d (%05d) "
                    "%02d -> %02d (%05d) "
                    "%02d <- %02d (%05d) "
                     "\n", 
              ThisTask, PrevTask, head_sendcount,
              ThisTask, PrevTask, head_recvcount,
              ThisTask, NextTask, tail_sendcount,
              ThisTask, NextTask, tail_recvcount
              );
        }
    }
    if(do_head) MPI_Wait(&head_req, MPI_STATUS_IGNORE);
    if(do_tail) MPI_Wait(&tail_req, MPI_STATUS_IGNORE);


    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Type_free(&partype);
    g_slice_free(Node, before);
    g_slice_free(Node, behind);

    tree_free();
    tree_build();
    domain_mark_complete();
}
