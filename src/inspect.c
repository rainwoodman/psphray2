#include <glib.h>
#include <mpi.h>

#include "commonblock.h"
#include "fckey.h"
#include "par.h"
#include "tree.h"
#include "domain.h"
void inspect_par(int color) {
    ParIter iter1;
    ParIter iter2;
    for(par_t * i1 = par_iter_init(&iter1, D[color].psys),
        * i2 = (par_iter_init(&iter2, D[color].psys), par_iter_next(&iter2)); 
                i1 && i2; 
            i1 = par_iter_next(&iter1), i2 = par_iter_next(&iter2)) {
        if(!fckey_cmp(&i1->fckey, &i2->fckey) < 0) {
            g_warning("%02d par unordered %ld and %ld " 
             FCKEY_FMT ", " FCKEY_FMT, 
             ThisTask, par_iter_last_index(&iter1), par_iter_last_index(&iter2),
            FCKEY_PRINT(i1->fckey), FCKEY_PRINT(i2->fckey));
        }
    }
    TAKETURNS {
        par_t * first = par_index(D[color].psys, 0);
        par_t * last = par_index(D[color].psys, -1);
            g_print("D%04d local Par(%ld): ID = %ld - %ld, "
              "KEY = " FCKEY_FMT " - " FCKEY_FMT "\n", 
              D[color].index,
              par_get_length(D[color].psys), 
              first->id, last->id,
              FCKEY_PRINT(first->fckey), FCKEY_PRINT(last->fckey));
//              inspect_tree();
    }
}
void inspect_tree(int color) {
    TreeIter iter;
    intptr_t count = 0;
    g_print("D%04d tree dump", D[color].index);
    for(Node * node = tree_iter_init(&iter, D[color].treestore, NULL);
        node;
        node = tree_iter_next(&iter)) {
        #if 0
        if(node->type != NODE_TYPE_LEAF)
        #endif
        g_print(NODE_FMT "\n", 
                NODE_PRINT(node[0]));
        count ++;
    }
    g_print("total %ld\n", count);
}

void inspect_domain_table() {
    for(int i = 0; i < NDomain; i++) {
        g_print("D%04d %04d:%02d " FCKEY_FMT "-" FCKEY_FMT "\n",
            i, DT[i].HostTask, DT[i].Color, 
            FCKEY_PRINT(DT[i].first),
            FCKEY_PRINT(DT[i].end));
    }
}
