#include <glib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <signal.h>
#include <execinfo.h>
#include "mpiu.h"
#include "paramfile.inc"
#include "commonblock.h"

extern void init_power(void);
extern void init_disp(void);
extern void init_filter(void);
extern void disp(int ax);
extern void filter(int ax, char* fname, int i);

extern double PowerSpec(double);
void abort() {
    MPI_Abort(MPI_COMM_WORLD, 1);
    /* shut up the compiler */
    exit(1);
}
void print_handler(const gchar * string) {
    fputs(string, stdout);
    fflush(stdout);
}

static char ** args = NULL;

static GOptionEntry entries[] =
{
  { "verbose", 'v', 0, G_OPTION_ARG_NONE, &CB.F.VERBOSE, "Be verbose", NULL },
  { "index", 'I', 0, G_OPTION_ARG_NONE, &CB.F.INDEX, "write the region map, but do not do the FFT, the default is to write delta and disp.", NULL },
  { "power", 'p', 0, G_OPTION_ARG_NONE, &CB.F.POWERONLY, "just write out the power spectrum", NULL },
  { G_OPTION_REMAINING, 0, 0, G_OPTION_ARG_STRING_ARRAY, &args, "", NULL },
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
    dinteger(keyfile, "IC", "SphereMode", &CB.IC.SphereMode, 1); 
    _double(keyfile, "IC", "BoxSize", &CB.BoxSize);
    _double(keyfile, "IC", "a", &CB.a);
    _double(keyfile, "IC", "Sigma8", &CB.C.Sigma8);
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
static int cmp_level_t_nmesh(level_t * l1, level_t *  l2) {
    return (l1->Nmesh > l2->Nmesh) - (l2->Nmesh > l1->Nmesh);
}
static void paramfile_read(char * filename) {
    GKeyFile * keyfile = g_key_file_new();
    GError * error = NULL;
    if(!g_key_file_load_from_file(keyfile, filename, G_KEY_FILE_KEEP_COMMENTS, &error)) {
        g_error("check param file %s:%s", filename, error->message);
    }

    _string(keyfile, "IO", "datadir", &CB.datadir);
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
        char * arg = g_key_file_get_string(keyfile, "Regions", keys[i], &error);
        if(!arg) {
            g_error("region failed :%s", error->message);
        }
        char ** sp = g_strsplit_set(arg, ",;", -1);
        int length = g_strv_length(sp);
        if(length != 4 && length != 6) {
            g_error("region must be provided with 4 or 6 numbers: x;y;z; sx[;sy;sz]");
        }
        for(int d = 0; d < 3; d++) {
            CB.IC.R[i].center[d] = atof(sp[d]);
            CB.IC.R[i].size[d] = atof(sp[length == 6 ? d + 3 : 3]);
        }
        g_strfreev(sp);
        g_free(arg);
        g_message("using region %d: center (%g %g %g) size (%g %g %g)", i,
                CB.IC.R[i].center[0], CB.IC.R[i].center[1], CB.IC.R[i].center[2], 
                CB.IC.R[i].size[0], CB.IC.R[i].size[1], CB.IC.R[i].size[2]);

    }
    g_strfreev(keys);

    keys = g_key_file_get_keys(keyfile, "Levels", &nkeys, &error);
    if(keys == NULL) {
        g_error("empty regions list:%s", error->message);
    }
    CB.IC.Scale = -1;
    CB.IC.NmeshAfter = INT_MAX;
    CB.IC.NLevels = g_strv_length(keys);
    CB.IC.Levels = g_new0(level_t, CB.IC.NLevels);
    for(int i = 0; keys[i]; i++) {
        char * args = g_key_file_get_string(keyfile, "Levels", keys[i], &error);
        char ** sp = g_strsplit_set(args, ";,", -1);
        int length = g_strv_length(sp);
        if(length < 2) {
            g_error("a level must have 2+ entries, Nmesh scaling factor, and optionally dm ptype");
        }
        int Nmesh = atoi(sp[0]);
        CB.IC.Levels[i].Nmesh = Nmesh;
        CB.IC.Levels[i].Scale = atof(sp[1]);
        CB.IC.Levels[i].use_longrangecut = 0;
        CB.IC.Levels[i].use_shortrangecut = 0;
        if(length == 4) {
            CB.IC.Levels[i].DownSample = atoi(sp[3]);
        } else {
            CB.IC.Levels[i].DownSample = 1;
        }

        if(Nmesh == CB.IC.Nmesh) {
            CB.IC.Scale = atof(sp[1]);
            if(length == 4) {
                CB.IC.DownSample = atoi(sp[3]);
            } else {
                CB.IC.DownSample = 1;
            }
        } else {
            if(Nmesh < CB.IC.Nmesh && Nmesh > CB.IC.NmeshBefore) {
                CB.IC.NmeshBefore = Nmesh;
                if(length == 4) {
                    CB.IC.DownSampleBefore = atoi(sp[3]);
                } else {
                    CB.IC.DownSampleBefore = 1;
                }
            }
            if(Nmesh > CB.IC.Nmesh && Nmesh < CB.IC.NmeshAfter) {
                CB.IC.NmeshAfter = Nmesh;
                if(length == 4) {
                    CB.IC.DownSampleAfter = atoi(sp[3]);
                } else {
                    CB.IC.DownSampleAfter = 1;
                }
            }
        }
        g_strfreev(sp);
        g_free(args);
    }
    qsort(CB.IC.Levels, CB.IC.NLevels, sizeof(level_t), cmp_level_t_nmesh);
    for(int i = 0; i < CB.IC.NLevels - 1; i++) {
        if(CB.IC.Levels[i+1].DownSample > 1) {
            CB.IC.Levels[i].use_shortrangecut = 1;
        }
    }
    for(int i = 1; i < CB.IC.NLevels; i++) {
        if(CB.IC.Levels[i - 1].use_shortrangecut) {
           CB.IC.Levels[i].use_longrangecut = 1;
        }
    }
    if(CB.IC.Scale == -1) {
        g_error("Nmesh (%d) does not present in paramfile ", CB.IC.Nmesh);
    }
    if(CB.IC.Scale == 0.0 && CB.F.INDEX) {
        g_error("no index to be made for the base level mesh (whose Scale == 0.0)");
    }
    g_message("using scale %g", CB.IC.Scale);
    g_strfreev(keys);

    char * data = g_key_file_to_data(keyfile, NULL, NULL);
    char * usedfilename = g_strdup_printf("%s/paramfile-used", CB.datadir);
    if(!g_file_set_contents(usedfilename, data, -1, &error)) {
        g_critical("failed to save used params to %s:%s", usedfilename, error->message);
        g_error_free(error);
    }
    g_free(usedfilename);
    g_free(data);
    g_key_file_free(keyfile);
}

static int decwidth(int n) {
    g_assert(n > 0);
    char * tmp = g_strdup_printf("%d", n-1);
    int width = strlen(tmp);
    g_free(tmp);
    return width;
}
int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);
    g_log_set_default_handler(log_handler, NULL);
    g_set_print_handler(print_handler);
    mpiu_init();

    struct sigaction sa;

    ROOTONLY {
        GError * error = NULL;
        GOptionContext * context = g_option_context_new("paramfile Nmesh");
        g_option_context_add_main_entries(context, entries, NULL);
        g_option_context_set_summary(context, 
"create the displacement and delta from a given IC by Fourier transforming"
"a random realization of the power spectrum. Disp and Delta of the selected regions (3-D boxes) are dumped into files [0-3].%d, where %d is the process id. 0, 1, 2 are the 3-vector of the displacement in code units (of KPC/h, usually), and 3 is the density delta. With -I meta data about which dumped point is written index is the index of the point in the regional box, and regions(char) is the regions.");
        if(!g_option_context_parse(context, &argc, &argv, &error)) {
            g_print("Option parsing failed: %s", error->message);
            abort();
        }
        if(g_strv_length(args) != 2) {
            g_print(g_option_context_get_help(context, FALSE, NULL));
            abort();
        }
        CB.IC.Nmesh = g_ascii_strtoll(args[1], NULL, 10);
        if(CB.IC.Nmesh == 0) {
            g_error("must give Nmesh!");
        }
        /* this will make use of CB.IC.Nmesh, thus later */
        g_message("Reading param file %s", args[0]);
        paramfile_read(args[0]);
        g_option_context_free(context);
    }
    common_block_sync();
    init_power();
    ROOTONLY {
        char * fname = g_strdup_printf("%s/power-%d.txt", CB.datadir, CB.IC.Nmesh);
        double K0 = 2 * G_PI / CB.BoxSize;
        FILE * fp = fopen(fname, "w");
        for(int i = 0; i < CB.IC.Nmesh; i++) {
            fprintf(fp, "%g %g\n", i * K0, PowerSpec(i * K0));
        }
        fclose(fp);
        free(fname);
    }
    if(CB.F.POWERONLY) {
        MPI_Finalize();
        return 0;
    }
    init_disp();
    init_filter();

    ROOTONLY {
        if(!CB.F.INDEX) {
            char * fname = g_strdup_printf("%s/meta-%d", CB.datadir, CB.IC.Nmesh);
            dump_filter(fname);
            free(fname);
        }
    }

    char * blocks[] = {"region", "index", "dispx", "dispy", "dispz", "delta"};
    int axstart, axend;
    if(CB.F.INDEX) {
        axstart = -2;
        axend = 0;
    } else {
        axstart = 0;
        axend = 4;
    }
    int width = decwidth(NTask);
    int dswidth = decwidth(CB.IC.DownSample);

    for(int ax = axstart; ax < axend; ax ++) {
        /* -2 is the region mask and -1 is the index map, they do not need the displacment field */
        char * fname = NULL;
        if(ax >= 0) disp(ax);
        for(int i = 0; i < CB.IC.DownSample; i++) {
            ROOTONLY {
                fname = g_strdup_printf("%s/%s-%d.header", CB.datadir, blocks[ax + 2], CB.IC.Nmesh);
                FILE * fp = fopen(fname, "w");
                fprintf(fp, "NTask = %d\n", NTask);
                fprintf(fp, "DownSample = %d\n", CB.IC.DownSample);
                fclose(fp);
                g_free(fname);
            }
            if(CB.IC.DownSample == 1) {
                fname = g_strdup_printf("%s/%s-%d.%0*d", CB.datadir, blocks[ax + 2], CB.IC.Nmesh, width, ThisTask);
            } else {
                fname = g_strdup_printf("%s/%s-%d.%0*d-%0*d", CB.datadir, blocks[ax + 2], CB.IC.Nmesh, dswidth, i, width, ThisTask);
            }
            filter(ax, fname, i);
        }
        g_free(fname);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    free_disp();
    MPI_Finalize();
}
