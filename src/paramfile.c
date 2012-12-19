#include <glib.h>
#include <math.h>
#include <stdint.h>
#include "commonblock.h"
#include "paramfile.inc"

static void SECTION_COSMOLOGY() {
    ddouble("Cosmology", "h", &CB.C.h, 0.72);
    ddouble("Cosmology", "H", &CB.C.H, 0.1);
    ddouble("Cosmology", "G", &CB.C.G, 43007.1);
    ddouble("Cosmology", "C", &CB.C.C, 3e5);
    ddouble("Cosmology", "OmegaB", &CB.C.OmegaB, 0.044);
    ddouble("Cosmology", "OmegaM", &CB.C.OmegaM, 0.26);
    ddouble("Cosmology", "OmegaL", &CB.C.OmegaL, 0.74);
}
static void SECTION_UNIT() {
    double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;

    ddouble("Unit", "UnitLength_in_cm", &UnitLength_in_cm, 3.085678e21); // 1 kpc/h
    ddouble("Unit", "UnitMass_in_g", &UnitMass_in_g, 1.989e43);  //1e10 solar masses/h
    ddouble("Unit", "UnitVelocity_in_cm_per_s", &UnitVelocity_in_cm_per_s, 1e5); //1 km/sec
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
    _string("IO", "datadir", &CB.datadir);
    _string("IO", "inputbase", &CB.inputbase);
    _string("IO", "snapbase", &CB.snapbase);
    _integer("IO", "SnapNumMajorBegin", &CB.SnapNumMajorBegin);
    dinteger("IO", "SnapNumMajorEnd", &CB.SnapNumMajorEnd, CB.SnapNumMajorBegin + 1);
    dinteger("IO", "NumReader", &CB.NReader, 0);
    dinteger("IO", "IDByteSize", &CB.IDByteSize, 8);

}

static void SECTION_DOMAIN() {
    ddouble("Domain", "MemImbalanceTol", &CB.MemImbalanceTol, 0.05);
    ddouble("Domain", "MemAllocFactor", &CB.MemAllocFactor, 2.0);
    dinteger("Domain", "SubDomainsPerTask", &CB.SubDomainsPerTask, 4);
}

static void SECTION_TREE() {
    dinteger("Tree", "NodeSplitThresh", &CB.NodeSplitThresh, 64);
}

static void SECTION_SFREFF() {
    ddouble("StarFormation", "CritOverDensity", &CB.SFREFF.CritOverDensity, 57.7);
    ddouble("StarFormation", "CritPhysDensity", &CB.SFREFF.CritPhysDensity, 0);
    ddouble("StarFormation", "FactorSN", &CB.SFREFF.FactorSN, 0.1);
    ddouble("StarFormation", "FactorEVP", &CB.SFREFF.FactorEVP, 1000.0);
    ddouble("StarFormation", "TempSupernova", &CB.SFREFF.TempSupernova, 1e8);
    ddouble("StarFormation", "TempClouds", &CB.SFREFF.TempClouds, 1e3);
    ddouble("StarFormation", "MaxSfrTimescale", &CB.SFREFF.MaxSfrTimescale, 1.5);
}
void paramfile_read(char * filename) {
    PARAMFILE_INIT(filename);

    SECTION_COSMOLOGY();
    SECTION_UNIT();
    SECTION_IO();
    SECTION_DOMAIN();
    SECTION_TREE();
    SECTION_SFREFF();

    char * data = g_key_file_to_data(keyfile, NULL, NULL);
    char * usedfilename = g_strdup_printf("%s/paramfile-used", CB.datadir);
    if(!g_file_set_contents(usedfilename, data, -1, &error)) {
        g_critical("failed to save used params to %s:%s", usedfilename, error->message);
        g_error_free(error);
    }
    g_free(usedfilename);
    g_free(data);
    PARAMFILE_FINAL();
}
