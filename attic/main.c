#include <glib.h>
#include <stdint.h>
	
#include "par.h"

int main(int argc, char * argv[]) {
    g_mem_set_vtable(glib_mem_profiler_table);
    register_ptype(0, "DM", 0, TRUE);
    register_ptype(1, "GAS", 24, FALSE);
    PStore * pstore = pstore_new(4);
    int i;
    g_random_set_seed(0);
    int NPAR = 1000000;
    for(i = 0; i < NPAR; i++) {
        ipos_t ipos[3];
        ipos[0] = g_random_int_range(0, IPOS_LIMIT);
        ipos[1] = g_random_int_range(0, IPOS_LIMIT);
        ipos[2] = g_random_int_range(0, IPOS_LIMIT);
    //    g_message("inserting %d", i);
        pstore_insert(pstore, ipos, 0);
    }
    for(i = 0; i < NPAR; i++) {
        ipos_t ipos[3];
        ipos[0] = g_random_int_range(0, IPOS_LIMIT);
        ipos[1] = g_random_int_range(0, IPOS_LIMIT);
        ipos[2] = g_random_int_range(0, IPOS_LIMIT);
    //    g_message("inserting %d", i);
        pstore_insert(pstore, ipos, 1);
    }
    pstore_check(pstore);
    for(i = 0; i < NPAR; i++) {
        int ind = g_random_int_range(0, pstore->root->size);
        Par * p = pstore_get_nearby(pstore, ind);
        pstore_remove(pstore, p);
        g_slice_free(Par, p);
    }
    pstore_check(pstore);
    Par * p = pstore->root->first;
    g_message(PAR_FMT, PAR_PRINT(p[0])); 
    g_message("total length = %d", g_slist_length(p));
    pstore_remove(pstore, p);
    g_message("removed length = %d", g_slist_length(p));
    p = pstore->root->first;
    g_message("total length = %d", g_slist_length(p));
    g_message(PAR_FMT, PAR_PRINT(p[0])); 
    //g_mem_profile();
    //g_slice_debug_tree_statistics();
    return 0;
}
