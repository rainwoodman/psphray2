#include <glib.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpiu.h"
#include "paramfile.inc"
#include "commonblock.h"

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
  { "Nmesh", 'n', 0, G_OPTION_ARG_INT, &CB.IC.Nmesh, "Nmesh", NULL },
  { "scale", 's', 0, G_OPTION_ARG_DOUBLE, &CB.IC.Scale, "scale region sizes", NULL },
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

static void SECTION_COSMOLOGY(GKeyFile * keyfile) {
    ddouble(keyfile, "Cosmology", "h", &CB.C.h, 0.72);
    ddouble(keyfile, "Cosmology", "H", &CB.C.H, 0.1);
    ddouble(keyfile, "Cosmology", "G", &CB.C.G, 43007.1);
    ddouble(keyfile, "Cosmology", "C", &CB.C.C, 3e5);
    ddouble(keyfile, "Cosmology", "OmegaB", &CB.C.OmegaB, 0.044);
    ddouble(keyfile, "Cosmology", "OmegaM", &CB.C.OmegaM, 0.26);
    ddouble(keyfile, "Cosmology", "OmegaL", &CB.C.OmegaL, 0.74);
}
static void SECTION_IC(GKeyFile * keyfile) {
    _integer(keyfile, "IC", "Seed", &CB.IC.Seed);
    dinteger(keyfile, "IC", "WhichSpectrum", &CB.IC.WhichSpectrum, 1); /*1 for EH, 3 for Etsu, 2 is from file and broken*/
    _double(keyfile, "IC", "BoxSize", &CB.BoxSize);
    _double(keyfile, "IC", "a", &CB.a);
    ddouble(keyfile, "IC", "PrimordialIndex", &CB.IC.PrimordialIndex, 0.96);
    ddouble(keyfile, "IC", "ShapeGamma", &CB.IC.ShapeGamma, 0.201); /* only for Efstathiou specturm */
}
static void SECTION_UNIT(GKeyFile * keyfile) {
    double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;

    ddouble(keyfile, "Unit", "UnitLength_in_cm", &UnitLength_in_cm, 3.085678e21); // 1 kpc/h
    ddouble(keyfile, "Unit", "UnitMass_in_g", &UnitMass_in_g, 1.989e43);  //1e10 solar masses/h
    ddouble(keyfile, "Unit", "UnitVelocity_in_cm_per_s", &UnitVelocity_in_cm_per_s, 1e5); //1 km/sec
    CB.U.CM_h = 1 / (UnitLength_in_cm);
    CB.U.CM = CB.U.CM_h * CB.C.h;
    CB.U.SECOND_h = 1 / (UnitLength_in_cm / UnitVelocity_in_cm_per_s);
    CB.U.SECOND = CB.U.SECOND_h * CB.C.h;

    CB.U.MYEAR_h = 86400 * 365.2433 * 1e6 * CB.U.SECOND_h;
    CB.U.KPC_h = CB.U.CM_h * 3.085678e21;
    CB.U.GRAM_h = 1.0 / UnitMass_in_g;
    CB.U.GRAM = CB.U.GRAM_h * CB.C.h;
    CB.U.SOLARMASS = CB.U.GRAM * 1.989e43;
    CB.U.PROTONMASS = CB.U.GRAM * 1.6726e-24;
}
static void paramfile_read(char * filename) {
    GKeyFile * keyfile = g_key_file_new();
    GError * error = NULL;
    if(!g_key_file_load_from_file(keyfile, filename, G_KEY_FILE_KEEP_COMMENTS, &error)) {
        g_error("check param file %s:%s", filename, error->message);
    }

    SECTION_COSMOLOGY(keyfile);
    SECTION_UNIT(keyfile);
    SECTION_IC(keyfile);
    gsize nkeys = 0;
    char ** keys = g_key_file_get_keys(keyfile, "Regions", &nkeys, &error);
    if(keys == NULL) {
        g_error("empty regions list:%s", error->message);
    }
    CB.IC.NRegions = nkeys;
    CB.IC.R = g_new0(region_t, nkeys);

    for(int i = 0; keys[i]; i++) {
        size_t length;
        double * list = g_key_file_get_double_list(keyfile, "Regions", keys[i], &length, &error);
        if(!list) {
            g_error("region failed :%s", error->message);
        }
        if(length != 4 && length != 6) {
            g_error("region must be provided with 4 or 6 numbers: x;y;z; sx[;sy;sz]");
        }
        for(int d = 0; d < 3; d++) {
            CB.IC.R[i].center[d] = list[d];
            CB.IC.R[i].size[d] = list[length == 6 ? d + 3 : 3];
        }
        g_free(list);
        g_message("using region : center (%g %g %g) size (%g %g %g)", 
                CB.IC.R[i].center[0], CB.IC.R[i].center[1], CB.IC.R[i].center[2], 
                CB.IC.R[i].size[0], CB.IC.R[i].size[1], CB.IC.R[i].size[2]);

    }
    g_strfreev(keys);
    g_key_file_free(keyfile);
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
    common_block_sync();
    init_power();
    init_disp();
    for(int ax=-1; ax < 3; ax ++) {
        disp(ax);
    }
    free_disp();
}