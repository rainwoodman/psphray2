#include <glib.h>
#include <stdint.h>
#include "mpiu.h"
#include "commonblock.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <fftw3-mpi.h>

/* from disp.c */
extern double * Disp;
extern intptr_t Local_nx;
extern intptr_t Local_x_start;

intptr_t Local_x_end;
intptr_t Local_x_start_new;
intptr_t Local_x_end_new;
intptr_t Local_nx_new;

#define PR(i, j, k) \
    (Disp[(((i)-Local_x_start)*Nsample+(j))*(2*(Nsample/2+1))+(k)])
#define PRNew(i, j, k) \
    (DispNew[(((i)-Local_x_start_new)*Nsample_new+(j))*(2*(Nsample_new/2+1))+(k)])

static int Nmesh;
static int Nsample;

static void layout(int x0, int x1, int * out) {
    /* out[i] is the process i-th slab is stored */
    int size = x1 - x0;
    int recvcount[NTask];
    int rdisp[NTask];
    MPI_Allgather(&size, 1, MPI_INT, 
            recvcount, 1, MPI_INT, MPI_COMM_WORLD);

    rdisp[0] = 0;
    for(int i = 1; i < NTask; i++) {
        rdisp[i] = rdisp[i - 1] + recvcount[i];
    }
    int buf[size];
    for(int i = 0; i < size; i++) {
        buf[i] = ThisTask;
    }
    MPI_Allgatherv(buf, size, MPI_INT,
            out, recvcount, rdisp, MPI_INT, MPI_COMM_WORLD);
}

void degrade(int Level, int newlevel) {
    Nmesh = L[Level].Nmesh;
    int factor = Nmesh / L[newlevel].Nmesh;
    int DownSample = L[Level].DownSample;
    g_assert(DownSample == 1);
    Nsample = Nmesh;
    g_assert(Nsample % factor == 0);
    int Nsample_new = Nsample / factor;

    size_t slab_size = 2 * (Nsample / 2 + 1L) * Nsample;
    size_t slab_size_new = 2 * (Nsample_new / 2 + 1L) * Nsample_new;

    Local_x_end = Local_x_start + Local_nx;
    Local_x_start_new = ThisTask * (Nsample_new) / NTask;
    Local_x_end_new = (ThisTask + 1)* (Nsample_new) / NTask;
    Local_nx_new = Local_x_end_new - Local_x_start_new;

    int OldLayout[Nsample];
    int NewLayout[Nsample_new];

    layout(Local_x_start, Local_x_end, OldLayout);
    layout(Local_x_start_new, Local_x_end_new, NewLayout);

    double * DispNew = fftw_malloc(sizeof(double) * 
            slab_size_new * Local_nx_new);

    memset(DispNew, 0, sizeof(double) * slab_size_new * Local_nx_new);

    /* this is the new slabs to be sent*/
    double * ExportSlabs[Nsample_new];
    memset(ExportSlabs, 0, sizeof(void*) * Nsample_new);
    for(int i = Local_x_start; i < Local_x_end; i++) {
        int inew = i / factor;
        if(inew < Local_x_start_new || inew >= Local_x_end_new) {
            /* export this*/ 
            if(ExportSlabs[inew] == NULL) {
                ExportSlabs[inew] = fftw_malloc(sizeof(double) * slab_size_new);
                memset(ExportSlabs[inew], 0, sizeof(double) * slab_size_new);
            }
            for(int j = 0; j < Nsample; j++) {
                int jnew = j / factor;
                for(int k = 0; k < Nsample; k++) {
                    int knew = k / factor;
                    ExportSlabs[inew]
                        [jnew * 2 * (Nsample_new / 2 + 1) + knew] 
                        += PR(i, j, k);
                }
            }
        } else {
            for(int j = 0; j < Nsample; j++) {
                int jnew = j / factor;
                for(int k = 0; k < Nsample; k++) {
                    int knew = k / factor;
                    PRNew(inew, jnew, knew) += PR(i, j, k);
                }
            }
        }
    }
    MPI_Request ExportRequest[Nsample_new];
    for(int i = 0; i < Nsample_new; i++) {
        if(ExportSlabs[i]) {
            MPI_Isend(ExportSlabs[i], slab_size_new, MPI_DOUBLE,
                    NewLayout[i], 99999, MPI_COMM_WORLD, &ExportRequest[i]);
        } else {
            ExportRequest[i] = MPI_REQUEST_NULL;
        }
    }
    int LastRecvFrom = -1;
    for(int i = Local_x_start_new; i < Local_x_end_new; i++) {
        for(int u = 0; u < factor; u ++) {
            int iold = i * factor + u;
            if(iold < Local_x_start || iold >= Local_x_end) {
                /* Import */ 

                /* do not recv twice */
                if(OldLayout[iold] == LastRecvFrom) continue;
                LastRecvFrom = OldLayout[iold];

                double * tmp = fftw_malloc(sizeof(double) * slab_size_new);
                MPI_Recv(tmp, slab_size_new, MPI_DOUBLE, 
                        OldLayout[iold], 99999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int j = 0; j < Nsample_new; j++) {
                    for(int k = 0; k < Nsample_new; k++) {
                        PRNew(i, j, k) += tmp[j * 2 * (Nsample_new / 2 + 1) + k];
                    }
                }
                fftw_free(tmp);
            }
        } 
    }
    for(intptr_t i = Local_x_start_new * slab_size_new; i < Local_x_end_new * slab_size_new; i ++) {
        DispNew[i] *= 1.0 / (factor * factor * factor);
    }
    fftw_free(Disp);
    Disp = DispNew;
    Local_x_start = Local_x_start_new;
    Local_nx = Local_nx_new;
}

