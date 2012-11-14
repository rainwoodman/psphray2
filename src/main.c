#include <glib.h>
#include <mpi.h>

#include "commonblock.h"


static char ** paramfilename = NULL;
static GOptionEntry entries[] =
{
  { "verbose", 'v', 0, G_OPTION_ARG_NONE, &CB.F.VERBOSE, "Be verbose", NULL },
  { G_OPTION_REMAINING, 0, 0, G_OPTION_ARG_FILENAME_ARRAY, &paramfilename, "Configratuion file", NULL },
  { NULL }
};

int main(int argc, char * argv[]) {
    GError * error = NULL;
    GOptionContext * context = g_option_context_new("paramfile");
    g_option_context_add_main_entries(context, entries, NULL);

    MPI_Init(&argc, &argv);
    common_block_bootstrap();

    ROOTONLY {
        if(!g_option_context_parse(context, &argc, &argv, &error)) {
            g_print("Option parsing failed: %s", error->message);
            abort(1);
        }
        if(g_strv_length(paramfilename) != 1) {
            g_print(g_option_context_get_help(context, FALSE, NULL));
            abort(1);
        }
        paramfile_read(paramfilename[0]);
        g_message("Reading param file %s", paramfilename[0]);
    }

    barrier();
    common_block_sync();
    barrier();

    g_message("MPI Task: %d of %d, datadir=%s", ThisTask, NTask, CB.datadir);
    barrier();

    for(CB.SnapNumMajor = CB.SnapNumMajorBegin;
        CB.SnapNumMajor < CB.SnapNumMajorEnd;
        CB.SnapNumMajor ++) {
        snapshot_prepare();
        snapshot_read();
        g_print("local Par(%d): ID = %ld - %ld\n", NPAR, PAR(0).id, PAR(-1).id);
    }
    barrier();
    MPI_Finalize();
    return 0;
}
