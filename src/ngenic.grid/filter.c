#include <glib.h>
#include "mpiu.h"
#include "commonblock.h"
#include <math.h>
#include <stdio.h>
/* from disp.c */
extern double * Disp;
extern intptr_t Local_nx;
extern intptr_t Local_x_start;

#define PR(i, j, k) Disp[(((i) - Local_x_start) * Nmesh + (j)) * (2 * (Nmesh / 2 + 1)) + (k)]
static double R0;
static int (*IBottom)[3];
static int (*ITop)[3];
static int (*ISize)[3];
static int Nmesh;
double F_Omega(double a);

void init_filter() {
    Nmesh = CB.IC.Nmesh;
    R0 = CB.BoxSize / Nmesh;
    IBottom = (int (*) [3])g_new0(int, 3 * CB.IC.NRegions);
    ITop = (int (*) [3]) g_new0(int, 3 * CB.IC.NRegions);
    ISize = (int (*) [3]) g_new0(int, 3 * CB.IC.NRegions);
    for(int r = 0; r < CB.IC.NRegions; r++) {
        for(int d = 0; d < 3; d++) {
            /* make sure the center is in the box.*/
            double c = remainder(CB.IC.R[r].center[d], CB.BoxSize);
            while(c < 0) c += CB.BoxSize;
            IBottom[r][d] = floor((c - CB.IC.R[r].size[d] * 0.5 * CB.IC.Scale) / R0) - 1;
            ITop[r][d] = ceil((c + CB.IC.R[r].size[d] * 0.5 * CB.IC.Scale) / R0)+ 1;
            ISize[r][d] = ITop[r][d] - IBottom[r][d];
        }
        ROOTONLY g_message("Region %d, %dx%dx%d", r, 
                ISize[r][0],
                ISize[r][1],
                ISize[r][2]);
    }
}

intptr_t filter0(int i, int j, int k, int r) {
    int ipos[] = {i, j, k};
    int inside = 0;
    int irelpos[3];
    for(int d = 0; d < 3; d++) {
        int dx = ipos[d] - ITop[r][d];
        while (dx < 0) dx += Nmesh;
        while (dx >= Nmesh) dx -= Nmesh;
        if(dx >= ISize[r][d]) return -1;
        irelpos[d] = dx;
    }
    return (((intptr_t)irelpos[0] * ISize[r][1]) + irelpos[1]) * ISize[r][2] + irelpos[2];
}

void filter(int ax, char * fname) {
    FILE * fp = fopen(fname, "w");
    const int BS = 1024 * 1024 * 8;
    char buffer[BS];
    char *be = &buffer[BS];
    char *bp = buffer;
    int i, j, k;
    for(i = 0; i < Nmesh; i ++) {
        if(i < Local_x_start) continue;
        if(i >= Local_x_start + Local_nx) continue;
        for(j = 0; j < Nmesh; j ++) {
            for(k = 0; k < Nmesh; k ++) {
                for(int r = 0; r < CB.IC.NRegions; r++) {
                    if(CB.IC.Scale == 0.0) {
                        float data = PR(i, j, k);
                        * (float * ) bp = data;
                        bp += sizeof(float);
                        break;
                    } else {
                        intptr_t index = filter0(i, j, k, r);
                        if(index >= 0) {
                            if(ax == -2) {
                                * (char*) bp = r;
                                bp ++;
                            } else if(ax == -1) {
                                * (intptr_t *) bp = index;
                                bp += sizeof(intptr_t);
                            } else {
                                float data = PR(i, j, k);
                                * (float * ) bp = data;
                                bp += sizeof(float);
                            }
                        }
                    }
                    if(bp == be) {
                        fwrite(buffer, BS, 1, fp);
                        bp = buffer;
                    }
                }
            }
        }
    }
    if(bp != buffer) {
        fwrite(buffer, bp - buffer, 1, fp);
    }
    fclose(fp);
}

void dump_filter(char * fname) {
    FILE * fp = fopen(fname, "w");
    double hubble_a = CB.C.H * sqrt(CB.C.OmegaM / pow(CB.a, 3) + 
            (1 - CB.C.OmegaM - CB.C.OmegaL) / pow(CB.a, 2) + CB.C.OmegaL);
    double vel_prefac = CB.a * hubble_a * F_Omega(CB.a);
    vel_prefac /= sqrt(CB.a);   /* converts to Gadget velocity */
    fprintf(fp, "vfact = %g\n", vel_prefac);
    fprintf(fp, "NTask = %d\n", NTask);
    fprintf(fp, "NRegion = %d\n", CB.IC.NRegions);
    for(int r = 0; r < CB.IC.NRegions; r++) {
        fprintf(fp, "Offset[%d] = [%d, %d, %d]\n", r, 
            IBottom[r][0], 
            IBottom[r][1], 
            IBottom[r][2]
        );
        fprintf(fp, "Size[%d] = [%d, %d, %d]\n", r, 
            ISize[r][0], 
            ISize[r][1], 
            ISize[r][2]
        );
    }
    fclose(fp);
}
