#include <glib.h>
#include <mpi.h>
#include <string.h>
#include "commonblock.h"
#include "fckey.h"
#include "par.h"
#include "tree.h"
#include "snapshot.h"
#include "domain.h"

Domain * D = NULL;
int NDomain = 0;
int NColor = 0;
DomainTable * DT = NULL;

static void cumbincount(fckey_t * POV, intptr_t * N) {
    for(int i = 0; i < NDomain - 1; i++) {
        N[i] = par_search_by_fckey(PAR_BUFFER_IN, &POV[i]);
    }
    N[NDomain - 1] = par_get_length(PAR_BUFFER_IN);
    MPI_Allreduce(MPI_IN_PLACE, N, NDomain, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
}
static void bincount_local(fckey_t * POV, intptr_t * N) {
    for(int i = 0; i < NDomain - 1; i++) {
        N[i] = par_search_by_fckey(PAR_BUFFER_IN, &POV[i]);
    }
    N[NDomain - 1] = par_get_length(PAR_BUFFER_IN);
    for(int i = NDomain - 1; i > 0; i--) {
        N[i] = N[i] - N[i - 1];
    }
}
/* POV will chop off the filling curve almost evenly */
static void find_pov(fckey_t * POV) {
    int i;
    fckey_t * POVleft =  g_new0(fckey_t, NDomain);
    fckey_t * POVright = g_new0(fckey_t, NDomain);
    fckey_t * POVmid = g_new0(fckey_t, NDomain);

    for(i = 0; i < NDomain; i++) {
        fckey_set_max(&POVright[i]);
        fckey_set_max(&POVmid[i]);
    }
    intptr_t * cumNtot_desired = g_new0(intptr_t, NDomain);
    intptr_t * cumNtot = g_new0(intptr_t, NDomain);

    intptr_t Ntot = 0;
    intptr_t NPARin = par_get_length(PAR_BUFFER_IN);

    MPI_Allreduce(&NPARin, &Ntot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    for(i = 0; i < NDomain; i++) {
        cumNtot_desired[i] = (i + 1)* Ntot / NDomain;
    }

    intptr_t tol = Ntot / NDomain * CB.MemImbalanceTol;
    int k = 0;
    while(-1) {
        for(i = 0; i < NDomain -1; i++) {
            fckey_center(&POVmid[i], &POVleft[i], &POVright[i]);
        }
        cumbincount(POVmid, cumNtot);
        for(i = 0; i < NDomain; i++) {
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
        for(i = 0; i < NDomain; i++) {
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
        for(i = 0; i < NDomain; i++) {
            g_debug("domain %d"
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
    for(i = 0; i < NDomain; i++) {
        POV[i] = POVmid[i];
    }

    g_free(cumNtot);
    g_free(cumNtot_desired);
    g_free(POVleft);
    g_free(POVright);
    g_free(POVmid);
}

/*
 * exchange items cross boundary, assuming
 * the initial ring topology of the domains
 * set up by domain_init.
 * */
static void exchange_cross_boundary(int elsize, void * first, void * last, void * before, void * behind) {
    MPI_Datatype item_type;
    MPI_Type_contiguous(elsize, MPI_BYTE, &item_type);
    MPI_Type_commit(&item_type);

    MPI_Barrier(MPI_COMM_WORLD);

    if(first && behind) {
        ROOTONLY {
            /* shift the first nodes before sending 
             * why? because  (assume two task and two colors)
             *
             * color  Task 0      Task 1
             *  0      0            1
             *  1      2            3
             *  2      4            5
             *   first(0,1)=2 -> 1=behind(1,0)
             *   first(0,2)=4 -> 3=behind(1,1)
             *   first(0,0)=0 -> 5=behind(1,2)
             * */
            void * temp = g_malloc0(elsize * NColor);
            memcpy(temp, (char*) first + elsize, elsize * (NColor - 1));
            memcpy((char*) temp + elsize * (NColor - 1), first, elsize);
            first = temp;
            /* free it later */
        }

        MPI_Sendrecv(first, NColor, item_type, PrevTask, 1,
            behind, NColor, item_type, NextTask, 1,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        ROOTONLY {
            /* because on root first array sent is a temporary storage
             * of the shifted */
            g_free(first);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(last && before) {
        MPI_Sendrecv(last, NColor, item_type, NextTask, 2,
            before, NColor, item_type, PrevTask, 2,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        ROOTONLY {
            /* shift the before nodes after receiving 
             * why? because  (assume two task and three colors)
             *
             * color  Task 0      Task 1
             *  0      0            1
             *  1      2            3
             *  2      4            5
             *   last(1, 2)=5 -> 0=before(0,0)
             *   last(1, 0)=1 -> 2=before(0,1)
             *   last(1, 1)=3 -> 4=before(0,2)
             * */
            void * temp = g_memdup(before, elsize * NColor);
            memcpy((char*) before + elsize, temp, elsize * (NColor - 1));
            memcpy(before, (char*) temp + elsize * (NColor - 1), elsize);
            g_free(temp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Type_free(&item_type);
}


void domain_init() {
    NColor = CB.SubDomainsPerTask;
    D = g_new0(Domain, NColor);
    NDomain = NColor * NTask;
    DT = g_new0(DomainTable, NColor * NTask);
    for(int i = 0; i < NTask; i++) {
        for(int color = 0; color < NColor; color++) {
            DT[color * NTask + i].HostTask = i;
            DT[color * NTask + i].Color = color;
            if(i == ThisTask) {
                D[color].index = color * NTask + i;
                D[color].prev = (color * NTask + i - 1 + NDomain) % NDomain;
                D[color].next = (color * NTask + i + 1) % NDomain;
                D[color].treestore = tree_store_alloc();
                D[color].psys = par_alloc();
                par_init(D[color].psys, "work");
            }
        }
    }
    #if 0
    /* Test exchange_cross_boundary */
    int * first = g_new0(int, NColor);
    int * last = g_new0(int, NColor);
    int * before = g_new0(int, NColor);
    int * behind = g_new0(int, NColor);
    for(int color = 0; color < NColor; color++) {
        first[color] = color;
        last[color] = color;
    }
    exchange_cross_boundary(sizeof(int), first, last, before, behind);
    for(int color = 0; color < NColor; color++) {
        g_print("color %d before %d, behind %d\n",
            color, before[color], behind[color]);
    }
    #endif
}
/*
 * decompose the domain->
 * 1) the input buffer is sorted. (will be freed after this)
 * 2) find POV points to cut the space filling curve, iteratively
 * 3) send particles to destinated Task's main buffers.
 * 4) sort the main buffer
 * 5) free the input buffer.
 * 6) adjust the domain boundary
 * 7) build the tree
 * */
void domain_decompose() {
    /* after domain_decompose, PAR is available */
    par_sort_by_fckey(PAR_BUFFER_IN);
    fckey_t * POV = g_new0(fckey_t, NDomain);
    intptr_t * N = g_new0(intptr_t, NDomain); /* particles per domain per task*/
    intptr_t * segN = g_new0(intptr_t, NDomain); /* particles per domain */

    find_pov(POV);

    /* segN will be used to allocate space for domains hosted locally */
    intptr_t Ntot = 0;
    bincount_local(POV, N);
    MPI_Allreduce(N, segN, NDomain, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    for(int i = 0; i < NDomain; i++) { Ntot += segN[i]; }


    /* report the load per task */
    ROOTONLY {
        /* no need to communicate as root has it all */
        intptr_t segNmax=-1, segNmin=-1;
        for(int i = 0; i < NTask; i++) {
            intptr_t seg = 0;
            for(int color = 0; color < NColor; color++) {
                seg += segN[i + color * NTask];
            }
            segNmax = MAX(segNmax, seg);
            segNmin = segNmin==-1?seg:MIN(segNmin, seg);
        }
        g_message("max load %ld, min load %ld", segNmax, segNmin);
    }

    /* now lets decide the communication layout */
    
    int *sendcounts = g_new0(int, NDomain);
    int *sdispls = g_new0(int, NDomain);

    /* how many from this Task will be sent to a domain,
     * and the offset */
    for(int i = 0; i < NDomain; i++) {
        sendcounts[i] = N[i];
        sdispls[i] = (i>0)?(sdispls[i-1] + sendcounts[i-1]):0;
    }

    MPI_Datatype partype = 0;
    MPI_Type_contiguous(sizeof(par_t), MPI_BYTE, &partype);
    MPI_Type_commit(&partype);

    /* sendcounts, recvcounts, sdispls rdispls are per domain, 
     * MPI communicate is done per Task. so we need to be clever.
     * 
     * The initial domain decomposition, each task contains SubDomainsPerTask
     * domains, and it is arranged in this way:
     * color 0: Task 0, 1, ....... NTask;
     * color 1: Task 0, 1, ....... NTask;
     * ......
     *
     * Thus we do MPI_Alltoall for each color.
     *
     * */
    for(int color = 0; color < NColor; color++) {
        int *recvcounts = g_new0(int, NTask);
        int *rdispls = g_new0(int, NTask);

        /* so we have a few slots before and after, 
         * to adjusted domain boundaries against the tree */
        intptr_t total = Ntot * CB.MemAllocFactor / NDomain;
        intptr_t before = (total - segN[D[color].index]) / 2;
        par_reserve(D[color].psys, total - before, before);
        par_t * recvbuffer = par_append(D[color].psys, segN[D[color].index]);

        MPI_Alltoall(&sendcounts[color * NTask], 1, MPI_INT,
                 recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

        for(int i = 0; i < NTask; i++) {
            rdispls[i] = (i>0)?(rdispls[i-1] + recvcounts[i-1]):0;
        }

        MPI_Alltoallv(par_index(PAR_BUFFER_IN, 0), &sendcounts[color * NTask], 
                    &sdispls[color * NTask], partype,
                    recvbuffer, recvcounts, 
                    rdispls, partype,
                    MPI_COMM_WORLD);


        g_free(rdispls);
        g_free(recvcounts);
    }
    MPI_Type_free(&partype);
    g_free(sdispls);
    g_free(sendcounts);
    g_free(segN);
    g_free(N);
    g_free(POV);
    par_destroy(PAR_BUFFER_IN);

    for(int color=0; color < NColor; color++) {
        par_sort_by_fckey(D[color].psys);
    }
}

static int node_intersects_domain(Node * node, DomainTable * d) {
    fckey_t tmp;
    tmp = node->key;
    fckey_set(&tmp, 3 * node->order);
#if 0
    g_print("xxxx " FCKEY_FMT ", " FCKEY_FMT " - ",
            FCKEY_PRINT(tmp), FCKEY_PRINT(*first));
    g_print("xxxx " FCKEY_FMT ", " FCKEY_FMT "\n",
            FCKEY_PRINT(node->key), FCKEY_PRINT(*last));
#endif
    /* inclusive */
    if(fckey_cmp(&tmp, &d->first) < 0) return 0;
    /* exculsive */
    if(fckey_cmp(&node->key, &d->end) >= 0) return 0;
    return 1;
}
static void update_ghosts() {
    fckey_t * first = g_new0(fckey_t, NColor);
    fckey_t * behind = g_new0(fckey_t, NColor);
    
    /* first will store the lower limit of fckeys 
     * hosted by the domain, and 
     * end will store the upper exclusive bound of fckeys,
     * so that as a whole, all domains cover the entire space */
    for(int color = 0; color < NColor; color++) {
        Node * root = tree_store_root(D[color].treestore);
        first[color] = tree_locate_down_fckey(D[color].treestore, root, 
                        &par_index(D[color].psys, 0)->fckey)->key;
        DT[D[color].index].first = first[color];
    }
    exchange_cross_boundary(sizeof(fckey_t), first, NULL, NULL, behind);
    for(int color = 0; color < NColor; color++) {
        DT[D[color].index].end = behind[color];
    } 
    g_free(first);
    g_free(behind);

    MPI_Datatype table_type;
    MPI_Type_contiguous(sizeof(DomainTable), MPI_BYTE, &table_type);
    MPI_Type_commit(&table_type);
    for(int color = 0; color < NColor; color++) {
        /* this is to avoid writing to the same location if the source
         * and dest are on the same hosting task, valgrind will discover this */
        DomainTable copy = DT[D[color].index];
        MPI_Allgather(&copy, 1, table_type,
            &DT[color * NTask], 1, table_type, MPI_COMM_WORLD);
    }
    MPI_Type_free(&table_type);
    /* first domain starts 0, last domain contains upto the tail */
    fckey_set_zero(&DT[0].first);
    fckey_set_max(&DT[NDomain-1].end);

    for(int color = 0; color < NColor; color++) {
        Node * leaf = NULL;
        intptr_t nleaf = 0;
        leaf = tree_store_get_leaf_nodes(D[color].treestore, &nleaf);
        
        for(intptr_t i = 0; i < nleaf; i++) {
            if(leaf[i].type != NODE_TYPE_GHOST) continue;
            Node * node = &leaf[i];
            int started = 0;
            for(int j = 0; j < NDomain; j++) {
                if(!started) {
                    if(node_intersects_domain(node, 
                        &DT[j])) {
                        node->ifirst = j;
                        started = 1;
                        node->npar = NDomain - j;
                        continue;
                    }
                } else {
                    if(!node_intersects_domain(node, 
                        &DT[j])) {
                        node->npar = j - node->ifirst;
                        break;
                    }
                }
            }
        }
    }

}
/*
 * Mark nodes that are completely within current Task.
 * node->complete is set to False if the another Task
 * overlaps this node.
 *
 * This is the last step of the domain decompostion.
 * called by domain_build_tree();
 * */

static void mark_complete() {
    fckey_t * first = g_new0(fckey_t, NColor);
    fckey_t * last = g_new0(fckey_t, NColor);
    fckey_t * before = g_new0(fckey_t, NColor);
    fckey_t * behind = g_new0(fckey_t, NColor);

    for(int color = 0; color < NColor; color++) {
        first[color] = par_index(D[color].psys, 0)->fckey;
        last[color] = par_index(D[color].psys, -1)->fckey;
    }
    exchange_cross_boundary(sizeof(fckey_t), first, last, before, behind);

    for(int color = 0; color < NColor; color++) {
        TreeIter iter;
        intptr_t complete_count = 0;
        intptr_t incomplete_count = 0;
        for(Node * node = tree_iter_init(&iter, D[color].treestore, NULL);
            node;
            node = tree_iter_next(&iter)) {
            if(node->type == NODE_TYPE_GHOST) {
                node->complete = FALSE;
                incomplete_count ++;
                continue;
            }
            if(D[color].index != 0 && 
                tree_node_contains_fckey(node, &before[color])) {
                node->complete = FALSE;
                incomplete_count ++;
              #if 1
                g_print("%04d incomplete" NODE_FMT ", key" FCKEY_FMT "\n",
                    D[color].index,
                    NODE_PRINT(node[0]), FCKEY_PRINT(before[0]));
              #endif
                continue;
            }
            if(D[color].index != NDomain -1 &&
                tree_node_contains_fckey(node, &behind[color])) {
                node->complete = FALSE;
                incomplete_count ++;
              #if 1
                g_print("%04d incomplete" NODE_FMT ", key" FCKEY_FMT "\n",
                    D[color].index,
                    NODE_PRINT(node[0]), FCKEY_PRINT(behind[0]));
              #endif
                continue;
            }
            node->complete = TRUE;
            complete_count ++;
        }
        TAKETURNS {
            g_print("%04d complete %ld incomplete %ld\n", 
                D[color].index,
                complete_count, incomplete_count);
        }
    }
    g_free(behind);
    g_free(before);
    g_free(first);
    g_free(last);
}

static void exchange_particles(Node * first, Node * last, Node * before, Node * behind) {
    MPI_Datatype partype;
    MPI_Type_contiguous(sizeof(par_t), MPI_BYTE, &partype);
    MPI_Type_commit(&partype);
    MPI_Request *frnt_req = g_new0(MPI_Request, NColor);
    MPI_Request *rear_req = g_new0(MPI_Request, NColor);

    intptr_t exchange = 0;
    for(int color=0; color < NColor; color++) {
        frnt_req[color] = MPI_REQUEST_NULL;
        rear_req[color] = MPI_REQUEST_NULL;

        int frnt_sendcount = 0;
        int frnt_recvcount = 0;
        int rear_sendcount = 0;
        int rear_recvcount = 0;

        int do_rear = tree_node_contains_node(&last[color], &behind[color]) || 
                    tree_node_contains_node(&behind[color], &last[color]);
        int do_frnt = tree_node_contains_node(&first[color], &before[color]) || 
                    tree_node_contains_node(&before[color], &first[color]);
        /* the domain with fewer particles send */
        if(do_rear) {
            if(last[color].npar < behind[color].npar) {
                /* send */
                rear_sendcount = last[color].npar;
                if(rear_sendcount) {
                    par_t * sendbuf = par_append(D[color].psys, -rear_sendcount);
                    MPI_Isend(sendbuf, rear_sendcount, partype,
                        DT[D[color].next].HostTask,
                        D[color].index,
                        MPI_COMM_WORLD, &rear_req[color]);
                }
            } else {
                /* recv */ 
                rear_recvcount = behind[color].npar;
                if(rear_recvcount) {
                    par_t * recvbuf = par_append(D[color].psys, rear_recvcount);
                    MPI_Irecv(recvbuf, rear_recvcount, partype, 
                        DT[D[color].next].HostTask,
                        D[color].next + NDomain,
                        MPI_COMM_WORLD, &rear_req[color]);
                }
            }
        }
        if(do_frnt) {
            if(first[color].npar > before[color].npar) {
                /* recv */
                frnt_recvcount = before[color].npar;
                if(frnt_recvcount) {
                    par_t * recvbuf = par_prepend(D[color].psys, frnt_recvcount);
                    MPI_Irecv(recvbuf, frnt_recvcount, partype, 
                        DT[D[color].prev].HostTask,
                        D[color].prev,
                        MPI_COMM_WORLD, &frnt_req[color]);
                }
            } else {
                /* send */
                frnt_sendcount = first[color].npar;
                if(frnt_sendcount) {
                    par_t * sendbuf = par_prepend(D[color].psys, -frnt_sendcount);
                    MPI_Isend(sendbuf, frnt_sendcount, partype,
                        DT[D[color].prev].HostTask,
                        D[color].index + NDomain,
                        MPI_COMM_WORLD, &frnt_req[color]);
                }
            }
        }
        exchange += frnt_sendcount + rear_sendcount + frnt_recvcount + rear_recvcount;
        TAKETURNS {
            if(frnt_sendcount + rear_sendcount + frnt_recvcount + rear_recvcount > 0) {
                g_print("%02d", ThisTask);
                if(frnt_sendcount)
                g_print(" %04df -> %04dr(% 5d)", D[color].index, D[color].prev, frnt_sendcount);
                if(frnt_recvcount)
                g_print(" %04df <- %04dr(% 5d)", D[color].index, D[color].prev, frnt_recvcount);
                if(rear_sendcount)
                g_print(" %04dr -> %04df(% 5d)", D[color].index, D[color].next, rear_sendcount);
                if(rear_recvcount)
                g_print(" %04dr <- %04df(% 5d)", D[color].index, D[color].next, rear_recvcount);
                g_print("\n");
            }
        }
    }
    
    /* now wait till the async calls finish */
    for(int color=0; color < NColor; color++) {
        MPI_Wait(&frnt_req[color], MPI_STATUS_IGNORE);
        MPI_Wait(&rear_req[color], MPI_STATUS_IGNORE);
    }

    intptr_t total_exchange = 0;

    MPI_Allreduce(&exchange, &total_exchange, 1, MPI_LONG, 
            MPI_SUM, MPI_COMM_WORLD);

    ROOTONLY {
        g_message("exchanged %ld particles ", total_exchange);
    }


    MPI_Type_free(&partype);
    g_free(frnt_req);
    g_free(rear_req);
}

/*
 * Finalize the domain decomposition and build a local tree
 * for the domain->
 *
 * We need to communicate once more to ensure the tree nodes
 * at the boundary of the domain are not chopped off.
 *
 * In the end, the tree nodes that are completely within the
 * local Task are marked with node->complete == True.
 *
 * We build the tree twice, first time is on the original
 * domain from domain_decompose; but this first tree is temporary,
 * we then decide which node to migrate to neighbour Tasks, migrate
 * the particles, and build the tree for the second time. 
 * The second tree is the real tree.
 * */

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

    Node * first = g_new0(Node, NColor);
    Node * last  = g_new0(Node, NColor);
    Node * before= g_new0(Node, NColor);
    Node * behind= g_new0(Node, NColor);

    /* first build a temporary tree */
    intptr_t skipped = 0;
    for(int color = 0; color < NColor; color++) {
        tree_store_init(D[color].treestore, D[color].psys, CB.NodeSplitThresh);
        skipped += tree_build(D[color].treestore);
    }
    MPI_Allreduce(MPI_IN_PLACE, &skipped, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    ROOTONLY {
        if(skipped > 0) {
            g_warning("particles out of box: %ld\n", skipped);
        }
    }

    for(int color = 0; color < NColor; color++) {
        Node * root = tree_store_root(D[color].treestore);
        first[color] = *tree_locate_down_fckey(D[color].treestore, root, &par_index(D[color].psys, 0)->fckey);
        last[color]  = *tree_locate_down_fckey(D[color].treestore, root, &par_index(D[color].psys,-1)->fckey);
    }
    /* exchange the boundary nodes of the trees */
    exchange_cross_boundary(sizeof(Node), first, last, before, behind);

    #if 0
    /* print the first last before exchange */
    for(int color = 0; color < NColor; color++) {
        TAKETURNS {
            g_print("%04d "
                " first " NODE_FMT 
                " last " NODE_FMT 
                " before " NODE_FMT 
                " behind " NODE_FMT "\n",
                D[color].index,
                NODE_PRINT(first[color]), 
                NODE_PRINT(last[color]), 
                NODE_PRINT(before[color]), 
                NODE_PRINT(behind[color]));
        }
    }
    #endif
    for(int color = 0; color < NColor; color++) {
        /* see if the node are joint, and exchange
         * the image nodes instead.
         * if each process is checking both first
         * and last, then we cover all four situations*/
        if(tree_node_contains_node(&behind[color], &last[color])) {
            last[color] = *tree_locate_up_image(D[color].treestore, &last[color], &behind[color]);
        }
        if(tree_node_contains_node(&before[color], &first[color])) {
            first[color] = *tree_locate_up_image(D[color].treestore, &first[color], &before[color]);
        }
    }

    exchange_cross_boundary(sizeof(Node), first, last, before, behind);

    /* after we have decided the nodes to exchange, the tree becomes irrelavant. 
     * the nodes that are useful are copied, and we only use the npar property. */
    for(int color=0; color < NColor; color++) {
        tree_store_destroy(D[color].treestore);
    }

    /* now we exchange the particles */
    exchange_particles(first, last, before, behind);

    g_free(before);
    g_free(behind);
    g_free(first);
    g_free(last);

    /* rebuild the tree */
    for(int color=0; color < NColor; color++) {
        tree_store_init(D[color].treestore, D[color].psys, CB.NodeSplitThresh);
        tree_build(D[color].treestore);
        /* create the ghost nodes ! mark_complete will fill
         * the right host task */
        tree_terminate(D[color].treestore);
    }
    mark_complete();
    update_ghosts();
}

void domain_cleanup() {
    for(int color = 0; color < NColor; color++) {
        par_destroy(D[color].psys);
        tree_store_destroy(D[color].treestore);
    }
}

void domain_destroy() {
    for(int color = 0; color < NColor; color++) {
        par_destroy(D[color].psys);
        par_free(D[color].psys);
        tree_store_destroy(D[color].treestore);
        tree_store_free(D[color].treestore);
    }
    g_free(D);
    D = NULL;
}
