#include <glib.h>
#include <gsl/gsl_rng.h>
#include <stdint.h>
#include "mpiu.h"
#include "commonblock.h"
#include <string.h>
#include <fftw3-mpi.h>
#include <math.h>

static fftw_complex * Cdata;
double * Disp;
intptr_t Local_nx;
intptr_t Local_x_start;
static fftw_plan InversePlan;
static size_t localsize;

static uint32_t * seedtable;
static int Nmesh, Nsample;
static double BoxSize;
static int SphereMode;
static double Dplus;
/* from powerspec.c */
extern double PowerSpec(double);
extern double GrowthFactor(double, double);

#define PK(i, j, k) Cdata[(((i) - Local_x_start) * Nmesh + (j)) * (Nmesh / 2 + 1) + (k)]

void init_disp() {
    Nmesh = CB.IC.Nmesh;
    Nsample = CB.IC.Nmesh;
    BoxSize = CB.BoxSize; 
    SphereMode = 1;

    Dplus = GrowthFactor(CB.a, 1.0);

    fftw_mpi_init();

    localsize = fftw_mpi_local_size_3d(Nmesh, Nmesh, Nmesh, MPI_COMM_WORLD, &Local_nx, &Local_x_start);

    gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    gsl_rng_set(random_generator, CB.IC.Seed);

    seedtable = g_new0(uint32_t, Nmesh * Nmesh);
    for(int i = 0; i < Nmesh / 2; i++)
    {
        int j;
        for(j = 0; j < i; j++)
            seedtable[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i + 1; j++)
            seedtable[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i; j++)
            seedtable[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i + 1; j++)
            seedtable[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i; j++)
            seedtable[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i + 1; j++)
            seedtable[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i; j++)
            seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i + 1; j++)
            seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }
    gsl_rng_free(random_generator);
}

void free_disp(void) {
    if(Cdata) free(Cdata);
    g_free(seedtable);
    fftw_destroy_plan(InversePlan);
    fftw_mpi_cleanup();
}

static int i_in_slab(int i) {
    return i >= Local_x_start && i < Local_x_start + Local_nx;
}

void fill_PK(int ax) {
    if(ax < 0) return;
    int i, ii, j, jj, k, kk;
    double K0 = 2 * G_PI / BoxSize;
    double fac = pow(K0, 1.5);
    double kvec[3];

    if(!Cdata) {
        /* allocate memory only if it is needed */
        Cdata = fftw_alloc_complex(localsize);
        if(!Cdata) {
            g_error("memory allocation failed");
        }
        ROOTONLY {
            g_message("allocated %td byte on Task %d for FFT's", localsize * sizeof(double) * 2, 0);
        }
        Disp = (double*) Cdata;
        /* In place b/c Cdata is an alias of Disp */
        InversePlan = fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, Cdata, Disp, MPI_COMM_WORLD, FFTW_ESTIMATE);
    }

    memset(Cdata, 0, Local_nx * Nmesh * (Nmesh / 2 + 1) * sizeof(fftw_complex));
    gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    for(i = 0; i < Nmesh; i++) {
        ii = (Nmesh - i) % Nmesh;
        if (!i_in_slab(i) & !i_in_slab(ii)) continue;

        for(j = 0; j < Nmesh; j++) {
            gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

            for(k = 0; k <= Nmesh / 2; k++) {
                double phase = gsl_rng_uniform(random_generator) * 2 * G_PI;
                double ampl;
                do {
                    ampl = gsl_rng_uniform(random_generator);
                } while(ampl == 0);

                /* singularity */
                if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2) continue;
                if(i == 0 && j == 0 && k == 0) continue;

                /* skip conjugates */
                if(k == 0) {
                    if(i == 0 && j > Nmesh / 2) continue; /* conjugate as 0, jj, 0 */
                    if(i > Nmesh / 2) continue;  /* conjugate as ii, jj, 0 */
                }

                if(i < Nmesh / 2) kvec[0] = i * K0;
                else kvec[0] = -(Nmesh - i) * K0;

                if(j < Nmesh / 2) kvec[1] = j * K0;
                else kvec[1] = -(Nmesh - j) * K0;

                /* the else branch never happens */
                if(k < Nmesh / 2) kvec[2] = k * K0;
                else kvec[2] = -(Nmesh - k) * K0;

                double kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
                double kmag = sqrt(kmag2);

                if(SphereMode == 1) {
                    if(kmag / K0 > Nsample / 2) continue;
                } else {
                    if(fabs(kvec[0]) / K0 > Nsample / 2)
                        continue;
                    if(fabs(kvec[1]) / K0 > Nsample / 2)
                        continue;
                    if(fabs(kvec[2]) / K0 > Nsample / 2)
                        continue;
                }

                double p_of_k = PowerSpec(kmag) * -log(ampl);
                double delta = fac * sqrt(p_of_k) / Dplus; /* scale back to starting redshift */

                double re, im;

                if(ax == 3) {
                    re = delta * cos(phase);
                    im = delta * sin(phase);
                } else if(ax >=0 && ax < 3) {
                    re = -kvec[ax] / kmag2 * delta * sin(phase);
                    im = kvec[ax] / kmag2 * delta * cos(phase);
                }

                if(k == 0) {
                    /* k=0 plane needs special treatment */
                    if(i == 0) {
                        /* PK(0, j, 0), PK(0, jj, 0) are conjugates */
                        jj = Nmesh - j; /* note: j!=0 surely holds at this point */

                        /* protect memory access */
                        i_in_slab(i)?(
                            PK(i, j, k)[0] = re,
                            PK(i, j, k)[1] = im,
                            PK(i, jj, k)[0] = re,
                            PK(i, jj, k)[1] = -im):0;
                    } else {

                        /* here comes i!=0 : conjugate can be on other processor! */
                        /* PK(i, j, 0), PK(ii, jj, 0) are conjugates */

                        jj = (Nmesh - j) % Nmesh;

                        /* protect memory access */
                        i_in_slab(i)?(
                            PK(i, j, k)[0] = re,
                            PK(i, j, k)[1] = im):0;
                        i_in_slab(ii)?(
                            PK(ii, jj, k)[0] = re,
                            PK(ii, jj, k)[1] = -im):0;
                    }
                } else {
                    /* k != 0, PK(i, j, ..) is independent */
                    i_in_slab(i)?(
                        PK(i, j, k)[0] = re,
                        PK(i, j, k)[1] = im):0;
                }
            /* end of k, j, i loops */
            }
        }
    }
    gsl_rng_free(random_generator);
}


void disp(int ax) {
    fill_PK(ax);
    fftw_execute(InversePlan); 
}
