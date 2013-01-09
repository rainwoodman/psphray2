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
intptr_t Local_ny;
intptr_t Local_y_start;
static fftw_plan InversePlan;
static size_t localsize;

static uint32_t * seedtable;
static int Nmesh, Nsample, DownSample, NmeshBefore;
static double BoxSize;
static int SphereMode;
static double Dplus;
/* from powerspec.c */
extern double PowerSpec(double);
extern double GrowthFactor(double, double);

#define PKT(i, j, k) Cdata[(((j) - Local_y_start) * Nsample + (i)) * (Nsample / 2 + 1) + (k)]

void init_disp() {
    Nmesh = CB.IC.Nmesh;
    NmeshBefore = CB.IC.NmeshBefore;
    DownSample = CB.IC.DownSample;
    Nsample = CB.IC.Nmesh / CB.IC.DownSample;
    BoxSize = CB.BoxSize; 
    SphereMode = 1;

    Dplus = GrowthFactor(CB.a, 1.0);

    fftw_mpi_init();

    localsize = fftw_mpi_local_size_3d_transposed(Nsample, Nsample,
             Nsample / 2+1, MPI_COMM_WORLD, 
            &Local_nx, &Local_x_start,
            &Local_ny, &Local_y_start);
}

void init_seedtable() {
    gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    gsl_rng_set(random_generator, CB.IC.Seed);

    seedtable = g_new0(uint32_t, Nsample * Nsample);
    uint32_t junk;
    #define ST(i, j) (((i) % DownSample == 0) && ((j) % DownSample == 0)?\
              &seedtable[(i) / DownSample * Nsample + (j) / DownSample]:\
              &junk)
    for(int i = 0; i < Nmesh / 2; i++)
    {
        int j;
        for(j = 0; j < i; j++) {
            ST(i, j)[0] = 0x7fffffff * gsl_rng_uniform(random_generator);
        }
        for(j = 0; j < i + 1; j++) {
            ST(j, i)[0] = 0x7fffffff * gsl_rng_uniform(random_generator);
        }
        for(j = 0; j < i; j++)
            ST((Nmesh - 1 - i), j)[0] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i + 1; j++)
            ST((Nmesh - 1 - j), i)[0] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i; j++)
            ST(i, (Nmesh - 1 - j))[0] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i + 1; j++)
            ST(j, (Nmesh - 1 - i))[0] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i; j++)
            ST((Nmesh - 1 - i), (Nmesh - 1 - j))[0] = 0x7fffffff * gsl_rng_uniform(random_generator);

        for(j = 0; j < i + 1; j++)
            ST((Nmesh - 1 - j), (Nmesh - 1 - i))[0] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }
    gsl_rng_free(random_generator);
}

void free_disp(void) {
    if(Cdata) free(Cdata);
    if(seedtable) g_free(seedtable);
    fftw_destroy_plan(InversePlan);
    fftw_mpi_cleanup();
}

static inline int j_in_slab(int j) {
    return j >= Local_y_start && j < Local_y_start + Local_ny;
}

void fill_PK(int ax) {
    if(ax < 0) return;
    int i, j, k;
    int dsi, dsii, dsj, dsjj, dsk, dskk;
    int kmod;
    double K0 = 2 * G_PI / BoxSize;
    double fac = pow(K0, 1.5);
    double kvec[3];

    if(!Cdata) {
        init_seedtable();
        /* allocate memory only if it is needed */
        Cdata = fftw_alloc_complex(localsize);
        if(!Cdata) {
            g_error("memory allocation failed");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        size_t llocalsize = localsize;
        MPI_Allreduce(MPI_IN_PLACE, &llocalsize, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
        ROOTONLY {
            g_message("max allocation %ld bytes for FFT, Nsample=%d", llocalsize * sizeof(double) * 2, Nsample);
        }
        Disp = (double*) Cdata;
        /* In place b/c Cdata is an alias of Disp */
        InversePlan = fftw_mpi_plan_dft_c2r_3d(Nsample, Nsample, Nsample, Cdata, Disp, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
        ROOTONLY {
            g_message("fft plan established");
        }
    }

    memset(Cdata, 0, localsize * sizeof(fftw_complex));

    gsl_rng * random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

    for(j = 0, dsj = 0; j < Nmesh; j += DownSample, dsj++) {
        dsjj = (Nsample - dsj) % Nsample;
        if (!j_in_slab(dsj) & !j_in_slab(dsjj)) continue;

        if (j == Nmesh / 2) continue;
        for(i = 0, dsi = 0; i < Nmesh; i += DownSample, dsi++) {
            if (i == Nmesh / 2) continue;
            dsii = (Nsample - dsi) % Nsample;

            gsl_rng_set(random_generator, seedtable[dsi * Nsample + dsj]);

            for(kmod = 0, dsk = 0, k = 0; 
                k <= Nmesh / 2; 
                k++, (++kmod==DownSample)?(kmod=0, dsk++):0) {
                double phase = gsl_rng_uniform(random_generator) * (2 * G_PI);
                double ampl;
                do {
                    ampl = gsl_rng_uniform(random_generator);
                } while(ampl == 0);

                /* to strictly preserve the random number sequence
                 * need to skip after phase and ampl are generated */

                /* singularity */
                if(k == Nmesh / 2) continue;
                if(i == 0 && j == 0 && k == 0) continue;
                /* not a sample point */
                if(kmod != 0) continue;

                /* skip conjugates */
                if(dsk == 0) {
                    /* conjugate of 0, jj, 0 */
                    if(dsi == 0 && dsj > Nsample / 2) continue; 
                    /* conjugate of ii, jj, 0 */
                    if(dsi > Nsample / 2) continue;  
                }

                /* use i, j, k for the K vector,
                 * and dsi, dsj, dsk for the memory location
                 * effectively this is doing a FFT on BoxSize/DownSample */
                if(i < Nmesh / 2) kvec[0] = i * K0;
                else kvec[0] = -(Nmesh - i) * K0;

                if(j < Nmesh / 2) kvec[1] = j * K0;
                else kvec[1] = -(Nmesh - j) * K0;

                /* the else branch never happens */
                if(k < Nmesh / 2) kvec[2] = k * K0;
                else kvec[2] = -(Nmesh - k) * K0;

                double kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
                double kmag = sqrt(kmag2);

                /* if this level has been downsampled, 
                 * remove the power from lower freqs, as they will
                 * be take from the lower level mesh.
                 * and assembled in convert.py. this has to be done
                 * to preserve the powerspectrum 
                 *
                 * In theory, for the no dowm sampling modes we shall
                 * also do this on the DownSampling = 1 meshes.
                 * however, the without a fourier interpolation(we do 
                 * nearest neighbour in conver.py)
                 * the density fluctuation won't be exactly the same;
                 * there is no need to pay this price.
                 * */
                if(SphereMode == 1) {
                    if(kmag / K0 > Nmesh / 2 ||
                    (DownSample > 1 && kmag / K0 >= NmeshBefore / 2)) 
                        continue;
                } else {
                    int d;
                    for(d = 0; d < 3; d++) {
                        if(fabs(kvec[d]) / K0 > Nmesh / 2 ||
                        (DownSample > 1 && 
                           fabs(kvec[d]) / K0 >= NmeshBefore / 2)) {
                            break;
                        }
                    }
                    /* if any of the dimension is out of bound */
                    if (d < 3) continue;
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

                if(dsk == 0) {
                    /* k=0 plane needs special treatment */
                    if(dsi == 0) {
                        /* PK(0, j, 0), PK(0, jj, 0) are conjugates */

                        /* protect memory access */
                        j_in_slab(dsj)?(
                            PKT(dsi, dsj, dsk)[0] = re,
                            PKT(dsi, dsj, dsk)[1] = im):0;
                        j_in_slab(dsjj)?(
                            PKT(dsi, dsjj, dsk)[0] = re,
                            PKT(dsi, dsjj, dsk)[1] = -im):0;
                    } else {

                        /* here comes j!=0 : conjugate can be on other processor! */
                        /* PK(i, j, 0), PK(ii, jj, 0) are conjugates */

                        /* protect memory access */
                        j_in_slab(dsj)?(
                            PKT(dsi, dsj, dsk)[0] = re,
                            PKT(dsi, dsj, dsk)[1] = im):0;
                        j_in_slab(dsjj)?(
                            PKT(dsii, dsjj, dsk)[0] = re,
                            PKT(dsii, dsjj, dsk)[1] = -im):0;
                    }
                } else {
                    /* k != 0, PK(i, j, ..) is independent */
                    j_in_slab(dsj)?(
                        PKT(dsi, dsj, dsk)[0] = re,
                        PKT(dsi, dsj, dsk)[1] = im):0;
                }
            /* end of k, j, i loops */
            }
        }
    }
    gsl_rng_free(random_generator);
}


void disp(int ax) {
    ROOTONLY g_message("filling axis %d", ax);
    fill_PK(ax);
    MPI_Barrier(MPI_COMM_WORLD);
    ROOTONLY g_message("fft axis %d", ax);
    fftw_execute(InversePlan); 
    ROOTONLY g_message("done fft axis %d", ax);
}
