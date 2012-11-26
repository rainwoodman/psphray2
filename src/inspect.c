#include <glib.h>
#include <mpi.h>

#include "commonblock.h"
void inspect_par(Domain * D) {
    for(intptr_t i = 0; NPAR(D) && i < NPAR(D) - 1; i++) {
        if(!fckey_cmp(&PAR(D, i).fckey, &PAR(D, i + 1).fckey) < 0) {
            g_warning("%02d par unordered %ld and %ld " 
             FCKEY_FMT ", " FCKEY_FMT, 
             ThisTask, i, i+1,
            FCKEY_PRINT(PAR(D, i).fckey), FCKEY_PRINT(PAR(D, i+1).fckey));
        }
    }
    TAKETURNS {
            g_print("%02d local Par(%ld): ID = %ld - %ld, "
              "KEY = " FCKEY_FMT " - " FCKEY_FMT "\n", 
              ThisTask, 
              NPAR(D), PAR(D, 0).id, PAR(D, -1).id,
              FCKEY_PRINT(PAR(D, 0).fckey), FCKEY_PRINT(PAR(D, -1).fckey));
//              inspect_tree();
    }
}
void inspect_tree(Domain * D) {
    TreeIter iter = {D->tree, NULL};
    Node * node = tree_iter_next(&iter);
    intptr_t count = 0;
    g_print("%02d tree dump\n", ThisTask);
    while(node) {
        #if 0
        if(node->type != NODE_TYPE_LEAF)
        #endif
        g_print(NODE_FMT "\n", 
                NODE_PRINT(node[0]));
        node = tree_iter_next(&iter);
        count ++;
    }
    g_print("%02d total %ld\n", ThisTask, count);
}
