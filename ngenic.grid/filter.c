#include <glib.h>
#include <stdint.h>
#include "mpiu.h"
#include "commonblock.h"
#include <math.h>
#include <stdio.h>
/* from disp.c */
extern double * Disp;
extern intptr_t Local_nx;
extern intptr_t Local_x_start;

#define PR(i, j, k) Disp[(((i) - Local_x_start) * Nsample + (j)) * (2 * (Nsample / 2 + 1)) + (k)]
static double R0;
static int (*IBottom)[3];
static int (*ITop)[3];
static int (*ISize)[3];
static intptr_t (*Stride)[3];
static int Nmesh;
static int Nsample;
static int DownSample;
static int Base;
double F_Omega(double a);

typedef struct {
    int region;
    intptr_t index;
} cross_t;

int sweep(cross_t * dest, int pos, int d, cross_t * src, int length) {
    /* sweep through all regions in cross, and see if the
     * d-th axis crosses the plane of pos. 
     * if not, the region will be removed from the list
     * if so, the raveled index of the position in the region 
     * will be updated.
     *
     * returns the new length of cross array.
     *
     * before calling this the first time, the cross array
     * shall be inialized to contain all regions and index=0
     * for each.
     * Stride and IBottom, ISize also need to be filled.
     * */
    if(Base) return 1;
    int newlength = 0;
    for(int i = 0; i < length; i++) {
        int r = src[i].region;
        int dx = pos - IBottom[r][d];
        while(dx < 0) dx += Nmesh;
        while(dx >= Nmesh) dx -= Nmesh;
        if(dx >= ISize[r][d]) continue;
        dest[newlength] = src[i];
        dest[newlength].index += Stride[r][d] * dx;
        newlength++;
    }
    return newlength;
}


void init_filter() {
    Nmesh = CB.IC.Nmesh;
    DownSample = CB.IC.DownSample;
    Nsample = Nmesh / DownSample;
    R0 = CB.BoxSize / Nmesh;
    Base = CB.IC.Scale == 0.0;
    IBottom = (int (*) [3])g_new0(int, 3 * CB.IC.NRegions);
    ITop = (int (*) [3]) g_new0(int, 3 * CB.IC.NRegions);
    ISize = (int (*) [3]) g_new0(int, 3 * CB.IC.NRegions);
    Stride = (intptr_t (*) [3]) g_new0(intptr_t, 3 * CB.IC.NRegions);
    for(int r = 0; r < CB.IC.NRegions; r++) {
        for(int d = 0; d < 3; d++) {
            /* make sure the center is in the box.*/
            double c = remainder(CB.IC.R[r].center[d], CB.BoxSize);
            while(c < 0) c += CB.BoxSize;
            IBottom[r][d] = floor((c - CB.IC.R[r].size[d] * 0.5 * CB.IC.Scale) / R0) - 1;
            ITop[r][d] = ceil((c + CB.IC.R[r].size[d] * 0.5 * CB.IC.Scale) / R0)+ 1;
            ISize[r][d] = ITop[r][d] - IBottom[r][d];
        }
        Stride[r][0] = ISize[r][1] * ISize[r][2];
        Stride[r][1] = ISize[r][2];
        Stride[r][2] = 1;
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
        int dx = ipos[d] - IBottom[r][d];
        while (dx < 0) dx += Nmesh;
        while (dx >= Nmesh) dx -= Nmesh;
        if(dx >= ISize[r][d]) return -1;
        irelpos[d] = dx;
    }
    return (((intptr_t)irelpos[0] * ISize[r][1]) + irelpos[1]) * ISize[r][2] + irelpos[2];
}

void filter(int ax, char * fname, int xDownSample) {
    FILE * fp = fopen(fname, "w");
    cross_t * cross0 = g_new0(cross_t, CB.IC.NRegions);
    cross_t * crossx = g_new0(cross_t, CB.IC.NRegions);
    cross_t * crossy = g_new0(cross_t, CB.IC.NRegions);
    cross_t * crossz = g_new0(cross_t, CB.IC.NRegions);
    /*copy all regions */
    int crossx_length;
    int crossy_length;
    int crossz_length;
    for(int icr = 0; icr < CB.IC.NRegions; icr++) {
        cross0[icr].region = icr;
        cross0[icr].index = 0;
    }
    const int BS = 1024 * 1024 * 8;
    char *buffer = g_malloc(BS);
    char *be = &buffer[BS];
    char *bp = buffer;
    int i, j, k, dsi, dsj, dsk;
    for(dsi = 0; dsi < Nsample; dsi ++) {
        i = dsi + xDownSample * Nsample;
        if(dsi < Local_x_start) continue;
        if(dsi >= Local_x_start + Local_nx) continue;
        crossx_length = sweep(crossx, i, 0, cross0, CB.IC.NRegions);
        if(!crossx_length) continue;
        for(j = 0; j < Nmesh; j ++) {
            crossy_length = sweep(crossy, j, 1, crossx, crossx_length);
            if(!crossy_length) continue;
            dsj = j % Nsample;
            for(k = 0; k < Nmesh; k ++) {
                crossz_length = sweep(crossz, k, 2, crossy, crossy_length);
                if(!crossz_length) continue;
                dsk = k % Nsample;
                if(Base) {
                    /* base scale, write displacement and delta,
                     * no special handling of ax==-2 and -1 */
                    float data = PR(dsi, dsj, dsk);
                    * (float * ) bp = data;
                    bp += sizeof(float);
                    if(bp == be) {
                        fwrite(buffer, BS, 1, fp);
                        bp = buffer;
                    }
                } else {
                    for(int icr = 0; icr < crossz_length; icr++) {
                        int r = crossz[icr].region;
                        intptr_t index = crossz[icr].index;
                        if(ax == -2) {
                            * (int*) bp = r;
                            bp += sizeof(int);
                        } else if(ax == -1) {
                            * (intptr_t *) bp = index;
                            bp += sizeof(intptr_t);
                        } else {
                            float data = PR(dsi, dsj, dsk);
                            * (float * ) bp = data;
                            bp += sizeof(float);
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
    g_free(cross0);
    g_free(crossz);
    g_free(crossy);
    g_free(crossx);
    g_free(buffer);
}

void dump_filter(char * fname) {
    FILE * fp = fopen(fname, "w");
    double hubble_a = CB.C.H * sqrt(CB.C.OmegaM / pow(CB.a, 3) + 
            (1 - CB.C.OmegaM - CB.C.OmegaL) / pow(CB.a, 2) + CB.C.OmegaL);
    double vel_prefac = CB.a * hubble_a * F_Omega(CB.a);
    vel_prefac /= sqrt(CB.a);   /* converts to Gadget velocity */
    fprintf(fp, "vfact = %g\n", vel_prefac);
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
