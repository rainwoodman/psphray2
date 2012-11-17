#include <glib.h>
#include <mpi.h>

#include "commonblock.h"

void abort() {
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static char ** paramfilename = NULL;
static GOptionEntry entries[] =
{
  { "verbose", 'v', 0, G_OPTION_ARG_NONE, &CB.F.VERBOSE, "Be verbose", NULL },
  { G_OPTION_REMAINING, 0, 0, G_OPTION_ARG_FILENAME_ARRAY, &paramfilename, "Configratuion file", NULL },
  { NULL }
};

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

    GError * error = NULL;
    GOptionContext * context = g_option_context_new("paramfile");
    g_option_context_add_main_entries(context, entries, NULL);

    MPI_Init(&argc, &argv);
    g_log_set_default_handler(log_handler, NULL);
    common_block_bootstrap();

    ROOTONLY {
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
    }

    MPI_Barrier(MPI_COMM_WORLD);
    common_block_sync();
    MPI_Barrier(MPI_COMM_WORLD);

    g_debug("MPI Task: %d of %d, datadir=%s", ThisTask, NTask, CB.datadir);
    MPI_Barrier(MPI_COMM_WORLD);

    for(CB.SnapNumMajor = CB.SnapNumMajorBegin;
        CB.SnapNumMajor < CB.SnapNumMajorEnd;
        CB.SnapNumMajor ++) {
        snapshot_read();
        par_sort_by_fckey(PAR_BUFFER_IN);
        domain_decompose();
        par_sort_by_fckey(PAR_BUFFER_MAIN);
        par_free(PAR_BUFFER_IN);
        tree_build();

        domain_adjust();

        tree_free();
        tree_build();
        domain_mark_complete();
        inspect_par();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

static void inspect_par() {
    for(intptr_t i = 0; NPAR && i < NPAR - 1; i++) {
        if(fckey_cmp(&PAR(i).fckey, &PAR(i + 1).fckey) > 0) {
            g_warning("%02d par unordered %ld and %ld " 
             FCKEY_FMT ", " FCKEY_FMT, 
             ThisTask, i, i+1,
            FCKEY_PRINT(PAR(i).fckey), FCKEY_PRINT(PAR(i+1).fckey));
        }
    }
    TAKETURNS {
            g_print("%02d local Par(%ld): ID = %ld - %ld, "
              "KEY = " FCKEY_FMT " - " FCKEY_FMT "\n", 
              ThisTask, 
              NPAR, PAR(0).id, PAR(-1).id,
              FCKEY_PRINT(PAR(0).fckey), FCKEY_PRINT(PAR(-1).fckey));
              inspect_tree();
    }
}
void inspect_tree() {
    TreeIter * iter = tree_iter_new(TREEROOT);
    Node * node = tree_iter_next(iter);
    intptr_t count = 0;
    g_print("%02d tree dump\n", ThisTask);
    while(node) {
        #if 0
        if(node->type != NODE_TYPE_LEAF)
        #endif
        g_print(NODE_FMT "\n", 
                NODE_PRINT(node[0]));
        node = tree_iter_next(iter);
        count ++;
    }
    g_print("%02d total %ld\n", ThisTask, count);
    tree_iter_free(iter);
}
