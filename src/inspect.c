#include <glib.h>
#include <mpi.h>

#include "commonblock.h"
void inspect_par(int color) {
    for(intptr_t i = 0; NPAR(color) && i < NPAR(color) - 1; i++) {
        if(!fckey_cmp(&PAR(color, i).fckey, &PAR(color, i + 1).fckey) < 0) {
            g_warning("%02d par unordered %ld and %ld " 
             FCKEY_FMT ", " FCKEY_FMT, 
             ThisTask, i, i+1,
            FCKEY_PRINT(PAR(color, i).fckey), FCKEY_PRINT(PAR(color, i+1).fckey));
        }
    }
    TAKETURNS {
            g_print("%02d %02d local Par(%ld): ID = %ld - %ld, "
              "KEY = " FCKEY_FMT " - " FCKEY_FMT "\n", 
              ThisTask, color,
              NPAR(color), PAR(color, 0).id, PAR(color, -1).id,
              FCKEY_PRINT(PAR(color, 0).fckey), FCKEY_PRINT(PAR(color, -1).fckey));
//              inspect_tree();
    }
}
void inspect_tree(int color) {
    TreeIter iter = {D[color].tree, NULL};
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
