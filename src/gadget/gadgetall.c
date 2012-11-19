#include <glib.h>
#include <stdint.h>
#include <math.h>
#include "commonblock.h"
#include "gadgetall.h"
struct global_data_all_processes All;
void init_gadget() {
    All.UnitMass_in_g = 1.0 / CB.U.GRAM_h;
    All.UnitLength_in_cm = 1.0 / CB.U.CM_h;
    All.UnitTime_in_s = 1.0 / CB.U.SECOND_h;
    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
    All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
    All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

    All.MinGasTemp = 5;

    InitCool(NULL);
    
    set_units_sfr();
    init_clouds();
}
