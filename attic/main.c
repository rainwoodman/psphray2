#include <glib.h>
#include <stdint.h>
#include <stdlib.h>
#include "mpiu.h"
#include "commonblock.h"
#include "par.h"
#include "hydro.h"
#include "snapshot.h"

static char ** paramfilename = NULL;
static GOptionEntry entries[] =
{
  { "verbose", 'v', 0, G_OPTION_ARG_NONE, &CB.F.VERBOSE, "Be verbose", NULL },
  { G_OPTION_REMAINING, 0, 0, G_OPTION_ARG_FILENAME_ARRAY, &paramfilename, "Configratuion file", NULL },
  { NULL }
};

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);
//    g_mem_set_vtable(glib_mem_profiler_table);
    mpiu_module_init();

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

    hydro_module_init();

    SnapHeader header;
    for(CB.SnapNumMajor = CB.SnapNumMajorBegin;
        CB.SnapNumMajor < CB.SnapNumMajorEnd;
        CB.SnapNumMajor ++) {
        PackedPar * pack = snapshot_read(&header);
        CB.a = header.a;
        ptrdiff_t i;
        uint64_t idmin = 0xfffffff, idmax = 0;
        ipos_t posmin[3] = {IPOS_LIMIT - 1, IPOS_LIMIT - 1, IPOS_LIMIT - 1} ,
               posmax[3] = {0, 0, 0};
        for(i = 0; i < pack->size; i++) {
            Par * par = pstore_pack_get(pack, i);
            if(par->id > idmax) idmax = par->id;
            if(par->id < idmin) idmin = par->id;
            if(ipos_compare(par->ipos, posmax) > 0) {
                posmax[0] = par->ipos[0];
                posmax[1] = par->ipos[1];
                posmax[2] = par->ipos[2];
            }
            if(ipos_compare(par->ipos, posmin) < 0) {
                posmin[0] = par->ipos[0];
                posmin[1] = par->ipos[1];
                posmin[2] = par->ipos[2];
            }
        }
        char * pmin = ipos_str(posmin), *pmax = ipos_str(posmax);
        g_message("%03d: Ntot = %td idmax = %lu, idmin = %lu"
                "posmin = %s posmax = %s",
                ThisTask, pack->size, idmax, idmin,pmin, pmax
                );
        domain(&pack);
        pstore_pack_free(pack);
    }
    g_mem_profile();
    //g_slice_debug_tree_statistics();
    MPI_Finalize();
    return 0;
}

