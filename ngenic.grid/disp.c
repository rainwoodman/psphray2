#include <glib.h>
#include <gsl/gsl_rng.h>
#include <stdint.h>
#include "mpiu.h"
#include "commonblock.h"
#include <string.h>
#include <fftw3-mpi.h>
#include <math.h>

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
static double Dplus;
/* from powerspec.c */
extern double PowerSpec(double);
extern double GrowthFactor(double, double);

#define Cdata ((fftw_complex *) Disp)
#define PKT(i, j, k) Cdata[(((j) - Local_y_start) * Nsample + (i)) * (Nsample / 2 + 1) + (k)]

void init_seedtable();
void init_disp(int Level) {
    Nmesh = L[Level].Nmesh;
    NmeshBefore = Level>0?L[Level-1].Nmesh:0;
    DownSample = L[Level].DownSample;
    Nsample = Nmesh / DownSample;
    BoxSize = CB.BoxSize; 
    Dplus = GrowthFactor(CB.a, 1.0);

    init_seedtable();


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
    if(Cdata) fftw_free(Cdata);
    fftw_destroy_plan(InversePlan);
}

static inline int j_in_slab(int j) {
    return j >= Local_y_start && j < Local_y_start + Local_ny;
}

void fill_PK(int ax) {
    if(ax < 0) return;
    int i, j, k;
    int dsi, dsii, dsj, dsjj, dsk, dskk;
    int kmod;
    /* What this does 
     * CFT = dx * iDFT (thus CFT has no 2pi factors and iCFT has, 
     *           same as wikipedia.)
     * iCFT(CFT) = (2pi) ** 3, thus iCFT also has no 2pi factors.
     * this is the same as the non-unitary transformation on wikipedia
     *
     * deltak = CFT(delta) = L**3 / N**3 iDFT(delta)
     * delta = L**-3 * DFT(deltak)
     * or,
     * delta = iCFT(deltak) * (2pi)**-3 = L**-3 DFT(deltak)
     *
     * Pk = CFT(xi) = deltak ** 2 * L ** -3
     * 
     * deltak = sqrt(Pk) * L ** 1.5
     *
     * delta = L **-1.5 DFT(sqrt(Pk * (2pi) **-3)) * (2pi) * 1.5
     *       = DFT(sqrt(Pk * (2pi) **-3) * (2pi /L) ** 1.5)
     *
     * observed:
     * delta = iDFT(sqrt(Pk * 2pi**-3) * (2pi/L) ** 1.5)
     * 
     * NOTE: PowerSpec is actually (Pk * 2pi ** -3).
     * DFT and iDFT differs by flipping sign of k. because input
     * is random this will just make it another realization.
     *
     * When DownSampling is enabled, effectively L is smaller.
     *
     * */
    double K0 = 2 * G_PI / BoxSize;
    double Kthresh = K0 * Nmesh / 2;
    double KthreshBefore = K0 * NmeshBefore / 2;
    /* Dplus scale back to starting redshift */
    /* DownSample ** 3 corrects for FFT<->CFT */
    double fac = pow(K0, 1.5) * pow(DownSample, 1.5) / Dplus; 

    double * ktable = g_new0(double, Nmesh);
    double * k2table = g_new0(double, Nmesh);
    for(i = 0; i < Nmesh; i ++) {
        if(i < Nmesh / 2) ktable[i] = i * K0;
        else ktable[i] = -(Nmesh - i) * K0;
        k2table[i] = ktable[i] * ktable[i];
    }

    /* allocate memory only if it is needed */
    Disp = (double*) fftw_alloc_complex(localsize);
    if(!Disp) {
        g_error("memory allocation failed");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    size_t llocalsize = localsize;
    MPI_Allreduce(MPI_IN_PLACE, &llocalsize, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    ROOTONLY {
        g_message("max allocation %ld bytes for FFT, Nsample=%d", llocalsize * sizeof(double) * 2, Nsample);
    }
    /* In place b/c Cdata is an alias of Disp */
    InversePlan = fftw_mpi_plan_dft_c2r_3d(Nsample, Nsample, Nsample, Cdata, Disp, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
    ROOTONLY {
        g_message("fft plan established");
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

                /* not a sample point */
                if(kmod != 0) continue;

                /* singularity */
                if(k == Nmesh / 2) continue;
                if(i == 0 && j == 0 && k == 0) continue;

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

                double kmag2 = k2table[i] + k2table[j] + k2table[k];
                double kmag = sqrt(kmag2);

                /* if this level has been downsampled, 
                 * remove the power from lower freqs, as they will
                 * be take from the lower level mesh.
                 * and assembled in convert.py. this has to be done
                 * to preserve the powerspectrum 
                 *
                 * In theory, for the no down sampling modes we shall
                 * also do this on the DownSampling = 1 meshes.
                 * however, the without a fourier interpolation(we do 
                 * nearest neighbour in conver.py)
                 * the density fluctuation won't be exactly the same;
                 * there is no need to pay this price.
                 * */
                if(CB.IC.SphereMode == 1) {
                    if(kmag > Kthresh ||
                            (DownSample > 1 && kmag <= KthreshBefore)) 
                        continue;
                } else {
                    if(ktable[i] > Kthresh ||
                            (DownSample > 1 && ktable[i] <= KthreshBefore))
                        continue;
                    if(ktable[j] > Kthresh ||
                            (DownSample > 1 && ktable[i] <= KthreshBefore))
                        continue;
                    if(ktable[k] > Kthresh ||
                            (DownSample > 1 && ktable[i] <= KthreshBefore))
                        continue;
                }

                /* PowerSpec is normalized to sigma8 ** 2 / (8 * pi**3) */
                double p_of_k = PowerSpec(kmag) * -log(ampl);
                double delta = fac * sqrt(p_of_k);
                double re, im;
                int dir[3] = {i, j, k};
                switch(ax) {
                    case 0:
                    case 1:
                    case 2:
                        re = sin(phase);
                        im = cos(phase);
                        delta *= ktable[dir[ax]] / kmag2;
                        re *= -delta;
                        im *= delta;
                        break;
                    case 3:
                        im = sin(phase);
                        re = cos(phase);
                        re *= delta;
                        im *= delta;
                        break;
                    default:
                        g_error("never reach here, ax is not 0 1 2 3");
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
    g_free(ktable);
    g_free(k2table);
}


void disp(int ax) {
    ROOTONLY g_message("filling axis %d", ax);
    fill_PK(ax);
    MPI_Barrier(MPI_COMM_WORLD);
    ROOTONLY g_message("fft axis %d", ax);
    fftw_execute(InversePlan); 
    ROOTONLY g_message("done fft axis %d", ax);
}

