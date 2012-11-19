#include <glib.h>
#include <math.h>
#include <stdint.h>
#include "commonblock.h"
#include "gadgetall.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
/*
 * This routine does cooling and star formation for
 * the effective multi-phase model.
 */

static double PhysDensThresh;
static double OverDensThresh;
static double EgySpecCold;
static double EgySpecSN;
double get_cloud_fraction(double density)
{
    double a3inv;
    double tsfr;
    double factorEVP, egyhot,  tcool, y, x;

    LoadIonizeParams(CB.a);
    a3inv = 1.0 / (CB.a * CB.a * CB.a);

    if(density * a3inv < PhysDensThresh) return 0;

    if(density < OverDensThresh) return 0;

    tsfr = sqrt(PhysDensThresh / (density * a3inv)) * CB.SFREFF.MaxSfrTimescale;

    factorEVP = pow(density * a3inv / PhysDensThresh, -0.8) * CB.SFREFF.FactorEVP;

    egyhot = EgySpecSN / (1 + factorEVP) + EgySpecCold;

    double ne = 1.0;
    tcool = GetCoolingTime(egyhot, density * a3inv, &ne);

    y = tsfr / tcool * egyhot / (CB.SFREFF.FactorSN * EgySpecSN - (1 - CB.SFREFF.FactorSN) * EgySpecCold);

    x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

    return x;
}

void init_clouds(void)
{
    double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;

    /* to be guaranteed to get z=0 rate */
    LoadIonizeParams(1.0);

    if(PhysDensThresh == 0)
    {
        A0 = CB.SFREFF.FactorEVP;

        egyhot = EgySpecSN / A0;

        meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

        u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
        u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


        dens = 1.0e6 * 3 * CB.C.H * CB.C.H / (8 * M_PI * CB.C.G);

        ne = 1.0;
        SetZeroIonization();
        tcool = GetCoolingTime(egyhot, dens, &ne);

        coolrate = egyhot / tcool / dens;

        x = (egyhot - u4) / (egyhot - EgySpecCold);

        PhysDensThresh =
            x / pow(1 - x,
                    2) * (CB.SFREFF.FactorSN * EgySpecSN - (1 -
                            CB.SFREFF.FactorSN) * EgySpecCold) /
                        (CB.SFREFF.MaxSfrTimescale * coolrate);

        if(ThisTask == 0)
        {
            g_message("A0= %g", A0);
            g_message("Computed: PhysDensThresh= %g (int units) %g h^2 cm^-3", 
                PhysDensThresh,
                PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
            g_message("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g", x);
            g_message("tcool=%g dens=%g egyhot=%g", tcool, dens, egyhot);
        }

    }
}


void set_units_sfr(void)
{
    double meanweight;

    OverDensThresh =
        CB.SFREFF.CritOverDensity * CB.C.OmegaB* 3 * CB.C.H * CB.C.H / (8 * M_PI * CB.C.G);
    PhysDensThresh = CB.SFREFF.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

    meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

    EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * CB.SFREFF.TempClouds;
    EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

    meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

    EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * CB.SFREFF.TempSupernova;
    EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

}



