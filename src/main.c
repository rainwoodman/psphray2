#include <glib.h>
#include <stdio.h>
#include "mpiu.h"

#include "commonblock.h"
#include "fckey.h"
#include "par.h"
#include "tree.h"
#include "snapshot.h"
#include "paramfile.h"
#include "domain.h"

void abort() {
    MPI_Abort(MPI_COMM_WORLD, 1);
    /* shut up the compiler */
    exit(1);
}
void print_handler(const gchar * string) {
    fputs(string, stdout);
    fflush(stdout);
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

    MPI_Init(&argc, &argv);
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

    init_gadget();
    g_message("MPI Task: %d of %d, datadir=%s", ThisTask, NTask, CB.datadir);
    MPI_Barrier(MPI_COMM_WORLD);

    domain_init();

    for(CB.SnapNumMajor = CB.SnapNumMajorBegin;
        CB.SnapNumMajor < CB.SnapNumMajorEnd;
        CB.SnapNumMajor ++) {
        SnapHeader h;
        snapshot_read(&h);
        CB.a = h.a;

        domain_decompose();
        domain_build_tree();
        for(int color=0; color < NColor; color++) {
            TAKETURNS {
                inspect_par(color);
                inspect_tree(color);
            }
        }
        ROOTONLY {
            inspect_domain_table();
        }
        domain_cleanup();
    }

    domain_destroy();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

