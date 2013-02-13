#ifndef __HYDRO_H__
#define __HYDRO_H__

#include "par.h"

typedef struct {
    real_t ie;
    real_t rho;
    real_t ye;
    real_t xHI;
    real_t sml;
    real_t sfr;
    real_t met;
} GasProp;
#define AS_GAS(par) ((GasProp *) ((Par *) (par))->prop)

typedef struct {
    real_t bhmass;
    real_t bhmdot;
} BHProp;
#define AS_BH(par) ((BHProp *) ((Par *) (par))->prop)

typedef struct {
    real_t sft;
    real_t met;
} StarProp;
#define AS_STAR(par) ((StarProp *) ((Par *) (par))->prop)

static void hydro_module_init() {
    register_ptype(0, "gas", sizeof(GasProp), FALSE);
    register_ptype(1, "dm", 0, TRUE);
    register_ptype(4, "star", sizeof(StarProp), FALSE);
    register_ptype(5, "bh", sizeof(BHProp), FALSE);
}
#endif
