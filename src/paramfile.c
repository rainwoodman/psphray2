#include <glib.h>
#include <mpi.h>
#include <math.h>
#include "commonblock.h"
static GError * error = NULL;
static GKeyFile * keyfile = NULL;

#define SCHEMA(name, type, mustval) \
static int _ ## name (char * group, char * key, type * value, type def) { \
    if(!g_key_file_has_group(keyfile, group) \
        || !g_key_file_has_key(keyfile, group, key, &error)) { \
        if(error != NULL) { \
            g_error("paramfile %s/%s:%s", group, key, error->message); \
        } \
        if(def == mustval) { \
            g_error("paramfile %s/%s is required", group, key); \
        } \
        g_key_file_set_ ## name(keyfile, group, key, def); \
        value[0] = def; \
        return FALSE; \
    } \
    value[0] = g_key_file_get_ ## name(keyfile, group, key, &error); \
    return TRUE; \
}
#define BADFLOAT g_strtod("NAN(BAD)",NULL)
#define BADINT 0xdeadbeef
SCHEMA(double, double, BADFLOAT);
SCHEMA(integer, int, BADINT);
SCHEMA(string, char *, NULL) ;

static void SECTION_COSMOLOGY() {
    _double("Cosmology", "h", &CB.C.h, 0.72);
    _double("Cosmology", "H", &CB.C.G, 0.1);
    _double("Cosmology", "G", &CB.C.G, 43007.1);
    _double("Cosmology", "C", &CB.C.C, 3e5);
    _double("Cosmology", "OmegaB", &CB.C.OmegaB, 0.044);
    _double("Cosmology", "OmegaM", &CB.C.OmegaM, 0.26);
    _double("Cosmology", "OmegaL", &CB.C.OmegaL, 0.74);
}
static void SECTION_UNIT() {
    double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;

    _double("Unit", "UnitLength_in_cm", &UnitLength_in_cm, 3.085678e21); // 1 kpc/h
    _double("Unit", "UnitMass_in_g", &UnitMass_in_g, 1.989e43);  //1e10 solar masses/h
    _double("Unit", "UnitVelocity_in_cm_per_s", &UnitVelocity_in_cm_per_s, 1e5); //1 km/sec
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

static void SECTION_IO() {
    _string("IO", "datadir", &CB.datadir, NULL);
    _string("IO", "inputbase", &CB.inputbase, NULL);
    _string("IO", "snapbase", &CB.snapbase, NULL);
    _integer("IO", "SnapNumMajorBegin", &CB.SnapNumMajorBegin, BADINT);
    _integer("IO", "SnapNumMajorEnd", &CB.SnapNumMajorEnd, BADINT);
    _integer("IO", "NumReader", &CB.NReader, BADINT);
    _integer("IO", "IDByteSize", &CB.IDByteSize, 8);

    if(CB.NReader > NTask) {
        g_warning("too many readers requested, using %d", NTask);
        CB.NReader = NTask;
    }
}

static void SECTION_DOMAIN() {
    _double("Domain", "MemImbalanceTol", &CB.MemImbalanceTol, 0.05);
    _double("Domain", "MemAllocFactor", &CB.MemAllocFactor, 2.0);
}

static void SECTION_TREE() {
    _integer("Tree", "NodeSplitThresh", &CB.NodeSplitThresh, 64);
}

void paramfile_read(char * filename) {
    keyfile = g_key_file_new();
    error = NULL;
    if(!g_key_file_load_from_file(keyfile, filename, 
        G_KEY_FILE_KEEP_COMMENTS, &error)) {
        g_error("check param file %s:%s", filename, error->message);
    }
    SECTION_COSMOLOGY();
    SECTION_UNIT();
    SECTION_IO();
    SECTION_DOMAIN();
    SECTION_TREE();

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
