#include <glib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpiu.h"
#include "commonblock.h"
#include "par.h"
static char ** paramfilename = NULL;
static GOptionEntry entries[] =
{
  { "verbose", 'v', 0, G_OPTION_ARG_NONE, &CB.F.VERBOSE, "Be verbose", NULL },
  { G_OPTION_REMAINING, 0, 0, G_OPTION_ARG_FILENAME_ARRAY, &paramfilename, "Configratuion file", NULL },
  { NULL }
};

#if 0
void abort() {
    MPI_Abort(MPI_COMM_WORLD, 1);
    /* shut up the compiler */
    exit(1);
}
#endif
void print_handler(const gchar * string) {
    fputs(string, stdout);
    fflush(stdout);
}

static void log_handler
    (const gchar *log_domain,
    GLogLevelFlags log_level,
    const gchar *message,
    gpointer unused_data) {

    g_log_default_handler(log_domain, log_level, message, unused_data);
    if((log_level & G_LOG_FLAG_FATAL) && NTask > 1) {
        g_on_error_stack_trace ("");
        abort();
    }
}

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    g_mem_set_vtable(glib_mem_profiler_table);
    g_log_set_default_handler(log_handler, NULL);
    g_set_print_handler(print_handler);
    mpiu_init();

    ROOTONLY {
        GError * error = NULL;
        GOptionContext * context = g_option_context_new("paramfile");
        g_option_context_add_main_entries(context, entries, NULL);

        if(!g_option_context_parse(context, &argc, &argv, &error)) {
            g_print("Option parsing failed: %s", error->message);
            abort();
        }
        if(g_strv_length(paramfilename) != 1) {
            g_print(g_option_context_get_help(context, FALSE, NULL));
            abort();
        }
        paramfile_read(paramfilename[0]);
        g_message("Reading param file %s", paramfilename[0]);
        g_option_context_free(context);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    common_block_sync();
    MPI_Barrier(MPI_COMM_WORLD);

    register_ptype(0, "DM", 0, TRUE);
    register_ptype(1, "GAS", 24, FALSE);

    PStore * pstore = pstore_new(4);
    int i;
    g_random_set_seed(0);
    int NPAR = 100000;
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
    PackedPar * pack = pstore_pack(pstore->root->first, pstore->root->size);
    pstore_pack_sort(pack);
    Par * par = pstore_unpack(pack);
    size_t c1 = 0;
    size_t c2 = 0;
    for(i = 0; i < pstore->root->size; i++) {
        if(par->type == 0) c1 ++;
        if(par->type == 1) c2 ++;
        if(par->next) {
            g_assert(ipos_compare(par->ipos, par->next->ipos) <= 0);
        }
        par = par->next;
    }
    g_assert(c1 == NPAR);
    g_assert(c2 == NPAR);
    pstore_check(pstore);
    for(i = 0; i < NPAR; i++) {
        int ind = g_random_int_range(0, pstore->root->size);
        Par * p = pstore_get_nearby(pstore, ind);
        pstore_remove(pstore, p);
        pstore_free_par(p);
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
