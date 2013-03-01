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
static int Nmesh;
static int Nsample;
static int DownSample;
static int Base;

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
        int dx = pos - R[r].ibottom[d];
        while(dx < 0) dx += Nmesh;
        while(dx >= Nmesh) dx -= Nmesh;
        if(dx >= R[r].isize[d]) continue;
        dest[newlength] = src[i];
        dest[newlength].index += R[r].stride[d] * dx;
        newlength++;
    }
    return newlength;
}


void init_filter() {
    Nmesh = CB.IC.Nmesh;
    DownSample = L[CB.IC.Level].DownSample;
    Nsample = Nmesh / DownSample;
    double Scale;
    Scale = L[CB.IC.Level].Scale;
    R0 = CB.BoxSize / Nmesh;
    Base = Scale == 0.0;
    for(int r = 0; r < NR; r++) {
        for(int d = 0; d < 3; d++) {
            /* make sure the center is in the box.*/
            double c = remainder(R[r].center[d], CB.BoxSize);
            while(c < 0) c += CB.BoxSize;
            R[r].ibottom[d] = floor((c - R[r].size[d] * 0.5 * Scale) / R0) - 1;
            R[r].itop[d] = ceil((c + R[r].size[d] * 0.5 * Scale) / R0)+ 1;
            R[r].isize[d] = R[r].itop[d] - R[r].ibottom[d];
        }
        R[r].stride[0] = R[r].isize[1] * R[r].isize[2];
        R[r].stride[1] = R[r].isize[2];
        R[r].stride[2] = 1;
        ROOTONLY g_message("Region %d, %dx%dx%d", r, 
                R[r].isize[0],
                R[r].isize[1],
                R[r].isize[2]);
    }

}

intptr_t filter0(int i, int j, int k, int r) {
    int ipos[] = {i, j, k};
    int inside = 0;
    int irelpos[3];
    for(int d = 0; d < 3; d++) {
        int dx = ipos[d] - R[r].ibottom[d];
        while (dx < 0) dx += Nmesh;
        while (dx >= Nmesh) dx -= Nmesh;
        if(dx >= R[r].isize[d]) return -1;
        irelpos[d] = dx;
    }
    return (((intptr_t)irelpos[0] * R[r].isize[1]) + irelpos[1]) * R[r].isize[2] + irelpos[2];
}

void filter(int ax, char * fname, int xDownSample) {
    FILE * fp = NULL;
    cross_t * cross0 = g_new0(cross_t, NR);
    cross_t * crossx = g_new0(cross_t, NR);
    cross_t * crossy = g_new0(cross_t, NR);
    cross_t * crossz = g_new0(cross_t, NR);
    /*copy all regions */
    int crossx_length;
    int crossy_length;
    int crossz_length;
    for(int icr = 0; icr < NR; icr++) {
        cross0[icr].region = icr;
        cross0[icr].index = 0;
    }
    const int BS = 1024 * 1024 * 8;
    char *buffer = g_malloc(BS);
    char *be = &buffer[BS];
    char *bp = buffer;
    int i, j, k, dsi, dsj, dsk;
    MPI_Barrier(MPI_COMM_WORLD);
    ROOTONLY g_message("filtering ax %d", ax);

    for(dsi = 0; dsi < Nsample; dsi ++) {
        i = dsi + xDownSample * Nsample;
        if(dsi < Local_x_start) continue;
        if(dsi >= Local_x_start + Local_nx) continue;
        crossx_length = sweep(crossx, i, 0, cross0, NR);
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
                        if(!fp) fp = fopen(fname, "w");
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
                        if(!fp) fp = fopen(fname, "w");
                        fwrite(buffer, BS, 1, fp);
                        bp = buffer;
                    }
                }
            }
        }
    }
    if(bp != buffer) {
        if(!fp) fp = fopen(fname, "w");
        fwrite(buffer, bp - buffer, 1, fp);
    }
    if(fp) fclose(fp);
    g_free(cross0);
    g_free(crossz);
    g_free(crossy);
    g_free(crossx);
    g_free(buffer);
    MPI_Barrier(MPI_COMM_WORLD);
    ROOTONLY g_message("done filtering ax %d", ax);
}

