#include <glib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <signal.h>
#include <execinfo.h>
#include "mpiu.h"
#include <string.h>
#include "paramfile.inc"
#include "commonblock.h"
#include <math.h>
#include <fftw3-mpi.h>
extern void init_power(void);
extern void init_filter(int Level);
extern void init_disp(int Level);
extern void init_seedtable(int Level);
extern void free_disp();
extern void disp(int ax);
extern void filter(int ax, char* fname, int i);
extern void degrade(int ax, int Level0, int Level1);
extern double F_Omega(double a);

extern double GrowthFactor(double, double);
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
  { "index", 'I', 0, G_OPTION_ARG_NONE, &CB.F.INDEX, "write the region map, but do not do the FFT, the default is to write index and delta and disp.", NULL },
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

int compare_logk(const void *a, const void *b) {
  if(((pow_table *) a)->logk < (((pow_table *) b)->logk))
    return -1;

  if(((pow_table *) a)->logk > (((pow_table *) b)->logk))
    return +1;

  return 0;
}


static void read_power_table() {
    FILE *fd;
    GList * ptr;
    double k, p;

    GList * table = NULL;

    g_message("using power from file %s", CB.IC.PowerSpectrumFile);
    if(!(fd = fopen(CB.IC.PowerSpectrumFile, "r"))) {
        g_error("can't read input spectrum in file '%s'\n", 
                CB.IC.PowerSpectrumFile);
    }

    while (fscanf(fd, " %lg %lg ", &k, &p) == 2) {
        struct pow_table * entry = g_new0(pow_table, 1);
        entry->logk = log10(k);
        entry->logD = log10(p);
        table = g_list_insert_sorted(table, (gpointer)entry, compare_logk);
    }
    fclose(fd);
    pow_table_alloc(g_list_length(table));
    int i;
    for(ptr = table, i = 0; ptr; ptr = ptr->next, i++) {
        memcpy(&PowerTable[i], ptr->data, sizeof(pow_table));
        printf("%g %g\n", PowerTable[i].logk, PowerTable[i].logD);
    }
    g_list_free_full(table, g_free);
    NPowerTable = i;
}


static void SECTION_COSMOLOGY(GKeyFile * keyfile) {
    ddouble(keyfile, "Cosmology", "h", &CB.C.h, 0.72);
    ddouble(keyfile, "Cosmology", "H", &CB.C.H, 0.1);
    ddouble(keyfile, "Cosmology", "G", &CB.C.G, 43007.1);
    ddouble(keyfile, "Cosmology", "C", &CB.C.C, 3e5);
    ddouble(keyfile, "Cosmology", "OmegaB", &CB.C.OmegaB, 0.044);
    ddouble(keyfile, "Cosmology", "OmegaM", &CB.C.OmegaM, 0.26);
    ddouble(keyfile, "Cosmology", "OmegaL", &CB.C.OmegaL, 0.74);
    _double(keyfile, "Cosmology", "Sigma8", &CB.C.Sigma8);
    dinteger(keyfile, "Cosmology", "WhichSpectrum", &CB.IC.WhichSpectrum, 1); /*1 for EH, 3 for Etsu, 2 is from file*/
    dstring(keyfile, "Cosmology", "PowerSpectrumFile", &CB.IC.PowerSpectrumFile, NULL); /*1 for EH, 3 for Etsu, 2 is from file*/
    ddouble(keyfile, "Cosmology", "PrimordialIndex", &CB.IC.PrimordialIndex, 0.96);
    ddouble(keyfile, "Cosmology", "ShapeGamma", &CB.IC.ShapeGamma, 0.201); /* only for Efstathiou specturm */
}
static void SECTION_IC(GKeyFile * keyfile) {
    _integer(keyfile, "IC", "Seed", &CB.IC.Seed);
    /* the primary mesh is where the fft is done.
     * coarse meshes are obtained by rebinning.*/
    _integer(keyfile, "IC", "NmeshPrimary", &CB.IC.NmeshPrimary);
    dinteger(keyfile, "IC", "SphereMode", &CB.IC.SphereMode, 1); 
    _double(keyfile, "IC", "BoxSize", &CB.BoxSize);
    _double(keyfile, "IC", "a", &CB.a);
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

    _string(keyfile, "IO", "datadir", &CB.datadir);
    SECTION_COSMOLOGY(keyfile);
    SECTION_UNIT(keyfile);
    SECTION_IC(keyfile);
    gsize nkeys = 0;

    char ** keys = g_key_file_get_keys(keyfile, "Regions", &nkeys, &error);
    if(keys == NULL) {
        g_error("empty regions list:%s", error->message);
    }
    g_assert(nkeys == 1);
    regions_alloc(nkeys);

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
            R[i].center[d] = atof(sp[d]);
            R[i].size[d] = atof(sp[length == 6 ? d + 3 : 3]);
        }
        g_strfreev(sp);
        g_free(arg);
        g_message("using region %d: center (%g %g %g) size (%g %g %g)", i,
                R[i].center[0], R[i].center[1], R[i].center[2], 
                R[i].size[0], R[i].size[1], R[i].size[2]);
    }

    g_strfreev(keys);
    keys = g_key_file_get_keys(keyfile, "Levels", &nkeys, &error);
    if(keys == NULL) {
        g_error("empty regions list:%s", error->message);
    }

    levels_alloc(g_strv_length(keys));

    for(int i = 0; keys[i]; i++) {
        char * args = g_key_file_get_string(keyfile, "Levels", keys[i], &error);
        char ** sp = g_strsplit_set(args, ";,", -1);
        int length = g_strv_length(sp);
        if(length < 2) {
            g_error("a level must have 2+ entries, Nmesh scaling factor, and optionally dm ptype");
        }
        int Nmesh = atoi(sp[0]);
        L[i].Nmesh = Nmesh;
        L[i].Scale = atof(sp[1]);
        if(length == 4) {
            L[i].DownSample = atoi(sp[3]);
        } else {
            L[i].DownSample = 1;
        }

        g_strfreev(sp);
        g_free(args);
    }

    levels_sort();

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

    if(CB.IC.WhichSpectrum == 2) {
        read_power_table();
    }
}

static int decwidth(int n) {
    g_assert(n > 0);
    char * tmp = g_strdup_printf("%d", n-1);
    int width = strlen(tmp);
    g_free(tmp);
    return width;
}

static void write_header(int Level, char * fname) {
    int DownSample = L[Level].DownSample;
    int Nmesh = L[Level].Nmesh;
    FILE * fp = fopen(fname, "w");
    fprintf(fp, "NTask = %d\n", NTask);
    fprintf(fp, "DownSample = %d\n", DownSample);
    double hubble_a = CB.C.H * sqrt(CB.C.OmegaM / pow(CB.a, 3) + 
        (1 - CB.C.OmegaM - CB.C.OmegaL) / pow(CB.a, 2) + CB.C.OmegaL);
    double vel_prefac = CB.a * hubble_a * F_Omega(CB.a);
    vel_prefac /= sqrt(CB.a);   /* converts to Gadget velocity */
    fprintf(fp, "vfact = %g\n", vel_prefac);
    fprintf(fp, "NRegion = %d\n", NR);
    fprintf(fp, "Dplus = %g\n", GrowthFactor(CB.a, 1.0));
    for(int r = 0; r < NR; r++) {
        fprintf(fp, "Offset[%d] = [%d, %d, %d]\n", r, 
            R[r].ibottom[0], 
            R[r].ibottom[1], 
            R[r].ibottom[2]
        );
        fprintf(fp, "Size[%d] = [%d, %d, %d]\n", r, 
            R[r].isize[0], 
            R[r].isize[1], 
            R[r].isize[2]
        );
    }
    fclose(fp);
}

void run_filter(int Level, int ax, int i) {
    static char * blocks[] = {"region", "index", "dispx", "dispy", "dispz", "delta"};
    int Nmesh = L[Level].Nmesh;
    int DownSample = L[Level].DownSample;

    int width = decwidth(NTask);
    int dswidth = decwidth(DownSample);
    char * fname;

    init_filter(Level); /* this will set R to the correct value*/

    ROOTONLY {
        fname = g_strdup_printf("%s/%s-%d.header", CB.datadir, blocks[ax + 2], Nmesh);
        write_header(Level, fname);
        free(fname);
    }
    if(DownSample == 1) {
        fname = g_strdup_printf("%s/%s-%d.%0*d", CB.datadir, blocks[ax + 2], Nmesh, width, ThisTask);
    } else {
        /* only power on this level is included, file postfix is '1' */
        fname = g_strdup_printf("%s/%s%s-%d.%0*d-%0*d", CB.datadir, blocks[ax + 2], (ax>=0)?"1":"", Nmesh, dswidth, i, width, ThisTask);
    }

    filter(ax, fname, i);

    g_free(fname);
}
int main(int argc, char * argv[]) {

    int Level;

    MPI_Init(&argc, &argv);
    g_log_set_default_handler(log_handler, NULL);
    g_set_print_handler(print_handler);

    mpiu_init();
    fftw_mpi_init();

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
        int Nmesh = g_ascii_strtoll(args[1], NULL, 10);
        if(Nmesh == 0) {
            g_error("must give Nmesh!");
        }
        /* this will make use of CB.IC.Nmesh, thus later */
        g_message("Reading param file %s", args[0]);
        paramfile_read(args[0]);
        g_option_context_free(context);

        if(Nmesh < CB.IC.NmeshPrimary) {
            g_error("Nmesh must be greater than NmeshPrimary in paramfile,"
                    "coarser meshes are produced by degrading the primary mesh.");
        }

        Level = levels_select(Nmesh);

        if(Level == 0 && CB.F.INDEX) {
         g_error("no index to be made for the base level mesh (whose Scale == 0.0)");
        }
        g_message("using scale %g", L[Level].Scale);
    }

    MPI_Bcast(&Level, sizeof(Level), MPI_BYTE, 0, MPI_COMM_WORLD);
    common_block_sync();

    int DownSample = L[Level].DownSample;
    int Nmesh = L[Level].Nmesh;

    init_power();

    ROOTONLY {
        char * fname = g_strdup_printf("%s/power-%d.txt", CB.datadir, Nmesh);
        double K0 = 2 * G_PI / CB.BoxSize;
        FILE * fp = fopen(fname, "w");
        for(int i = 0; i < Nmesh; i++) {
            fprintf(fp, "%g %g\n", i * K0, PowerSpec(i * K0));
        }
        fclose(fp);
        free(fname);
    }
    if(CB.F.POWERONLY) {
        MPI_Finalize();
        return 0;
    }

    int axstart, axend;
    if(CB.F.INDEX) {
        axstart = -2;
        axend = 0;
    } else {
        axstart = -2;
        axend = 4;
    }

    init_seedtable(Level);

    for(int ax = axstart; ax < axend; ax ++) {
        /* -2 is the region mask and -1 is the index map, they do not need the displacment field */
        init_disp(Level);
        disp(ax);
        if (Nmesh == CB.IC.NmeshPrimary) {
            /* this is the primary mesh; 
             * degrade to get the coarse meshes, too */
            for(int i = Level; i >= 0; i--) {
                MPI_Barrier(MPI_COMM_WORLD);
                ROOTONLY
                    g_message("running on level %d", L[i].Nmesh);
                if(i != Level)
                    degrade(ax, i + 1, i);
                /* these coarse levels never have DownSample,
                 * 0 or anything will work */
                run_filter(i, ax, 0);
            }
        } else {
            for(int i = 0; i < DownSample; i++) {
                run_filter(Level, ax, i);
            }
        }
        free_disp();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    fftw_mpi_cleanup();
    MPI_Finalize();
}
