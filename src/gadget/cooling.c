#include <glib.h>
#include <math.h>
#include "mpiu.h"
#include <stdio.h>
#include "commonblock.h"
#include "gadgetall.h"
#define NCOOLTAB  2000

#define SMALLNUM 1.0e-60
#define COOLLIM  0.1
#define HEATLIM	 20.0


static double Time;
static double XH = HYDROGEN_MASSFRAC;	/* hydrogen abundance by mass */
static double yhelium;

#define eV_to_K   11606.0
#define eV_to_erg 1.60184e-12


static double mhboltz;		/* hydrogen mass over Boltzmann constant */
static double ethmin;		/* minimum internal energy for neutral gas */

static double Tmin = 0.0;	/* in log10 */
static double Tmax = 9.0;
static double deltaT;

static double *BetaH0, *BetaHep, *Betaff;
static double *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp;
static double *GammaeH0, *GammaeHe0, *GammaeHep;

static double J_UV = 0, gJH0 = 0, gJHep = 0, gJHe0 = 0, epsH0 = 0, epsHep = 0, epsHe0 = 0;

static double ne, necgs, nHcgs;
static double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
static double gJH0ne, gJHe0ne, gJHepne;
static double nH0, nHp, nHep, nHe0, nHepp;



static double DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;


/* returns cooling time. 
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime(double u_old, double rho, double *ne_guess)
{
    double u;
    double ratefact;
    double LambdaNet, coolingtime;

    DoCool_u_old_input = u_old;
    DoCool_rho_input = rho;
    DoCool_ne_guess_input = *ne_guess;

    rho *= All.UnitDensity_in_cgs * CB.C.h * CB.C.h;	/* convert to physical cgs units */
    u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;


    nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho;

    u = u_old;

    LambdaNet = CoolingRateFromU(u, rho, ne_guess);

    /* bracketing */

    if(LambdaNet >= 0)		/* ups, we have actually heating due to UV background */
        return 0;

    coolingtime = u_old / (-ratefact * LambdaNet);

    coolingtime *= CB.C.h / All.UnitTime_in_s;

    return coolingtime;
}


void cool_test(void)
{
    double uin, rhoin, tempin, muin, nein;

    uin = 6.01329e+09;
    rhoin = 7.85767e-29;
    tempin = 34.0025;
    muin = 0.691955;

    nein = (1 + 4 * yhelium) / muin - (1 + yhelium);

    g_message("%g", convert_u_to_temp(uin, rhoin, &nein));
}


/* this function determines the electron fraction, and hence the mean 
 * molecular weight. With it arrives at a self-consistent temperature.
 * Element abundances and the rates for the emission are also computed
 */
double convert_u_to_temp(double u, double rho, double *ne_guess)
{
    double temp, temp_old, temp_new, max = 0, ne_old;
    double mu;
    int iter = 0;

    double u_input, rho_input, ne_input;

    u_input = u;
    rho_input = rho;
    ne_input = *ne_guess;

    mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
    temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

    do
    {
        ne_old = *ne_guess;

        find_abundances_and_rates(log10(temp), rho, ne_guess);
        temp_old = temp;

        mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);

        temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

        max =
            fmax(max,
                    temp_new / (1 + yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

        temp = temp_old + (temp_new - temp_old) / (1 + max);
        iter++;

        if(iter > (MAXITER - 10))
            g_warning("-> temp= %g ne=%g", temp, *ne_guess);
    }
    while(fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);

    if(iter >= MAXITER)
    {
        g_error("failed to converge in convert_u_to_temp()\n"
                "u_input= %g rho_input=%g ne_input=%g\n"
                "DoCool_u_old_input=%g DoCool_rho_input= %g \n"
                "DoCool_dt_input= %g DoCool_ne_guess_input= %g\n",
                u_input, rho_input, ne_input,
                DoCool_u_old_input, DoCool_rho_input, 
                DoCool_dt_input, DoCool_ne_guess_input);
    }

    return temp;
}


/* this function computes the actual abundance ratios 
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess)
{
    double neold, nenew;
    int j, niter;
    double Tlow, Thi, flow, fhi, t;

    double logT_input, rho_input, ne_input;

    logT_input = logT;
    rho_input = rho;
    ne_input = *ne_guess;

    if(logT <= Tmin)		/* everything neutral */
    {
        nH0 = 1.0;
        nHe0 = yhelium;
        nHp = 0;
        nHep = 0;
        nHepp = 0;
        ne = 0;
        *ne_guess = 0;
        return;
    }

    if(logT >= Tmax)		/* everything is ionized */
    {
        nH0 = 0;
        nHe0 = 0;
        nHp = 1.0;
        nHep = 0;
        nHepp = yhelium;
        ne = nHp + 2.0 * nHepp;
        *ne_guess = ne;		/* note: in units of the hydrogen number density */
        return;
    }

    t = (logT - Tmin) / deltaT;
    j = (int) t;
    Tlow = Tmin + deltaT * j;
    Thi = Tlow + deltaT;
    fhi = t - j;
    flow = 1 - fhi;

    if(*ne_guess == 0)
        *ne_guess = 1.0;

    nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

    ne = *ne_guess;
    neold = ne;
    niter = 0;
    necgs = ne * nHcgs;

    /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
    do
    {
        niter++;

        aHp = flow * AlphaHp[j] + fhi * AlphaHp[j + 1];
        aHep = flow * AlphaHep[j] + fhi * AlphaHep[j + 1];
        aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[j + 1];
        ad = flow * Alphad[j] + fhi * Alphad[j + 1];
        geH0 = flow * GammaeH0[j] + fhi * GammaeH0[j + 1];
        geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[j + 1];
        geHep = flow * GammaeHep[j] + fhi * GammaeHep[j + 1];

        if(necgs <= 1.e-25 || J_UV == 0)
        {
            gJH0ne = gJHe0ne = gJHepne = 0;
        }
        else
        {
            gJH0ne = gJH0 / necgs;
            gJHe0ne = gJHe0 / necgs;
            gJHepne = gJHep / necgs;
        }

        nH0 = aHp / (aHp + geH0 + gJH0ne);	/* eqn (33) */
        nHp = 1.0 - nH0;		/* eqn (34) */

        if((gJHe0ne + geHe0) <= SMALLNUM)	/* no ionization at all */
        {
            nHep = 0.0;
            nHepp = 0.0;
            nHe0 = yhelium;
        }
        else
        {
            nHep = yhelium / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
            nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
            nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
        }

        neold = ne;

        ne = nHp + nHep + 2 * nHepp;	/* eqn (38) */
        necgs = ne * nHcgs;

        if(J_UV == 0)
            break;

        nenew = 0.5 * (ne + neold);
        ne = nenew;
        necgs = ne * nHcgs;

        if(fabs(ne - neold) < 1.0e-4)
            break;

        if(niter > (MAXITER - 10))
            g_warning("ne= %g  niter=%d", ne, niter);
    }
    while(niter < MAXITER);

    if(niter >= MAXITER)
    {
        g_error("no convergence reached in find_abundances_and_rates()\n"
                "logT_input= %g  rho_input= %g  ne_input= %g\n"
                "DoCool_u_old_input=%g DoCool_rho_input= %g\n"
                "DoCool_dt_input= %g DoCool_ne_guess_input= %g\n",
                logT_input, rho_input, ne_input,
                DoCool_u_old_input, DoCool_rho_input, 
                DoCool_dt_input, DoCool_ne_guess_input);
    }

    bH0 = flow * BetaH0[j] + fhi * BetaH0[j + 1];
    bHep = flow * BetaHep[j] + fhi * BetaHep[j + 1];
    bff = flow * Betaff[j] + fhi * Betaff[j + 1];

    *ne_guess = ne;
}



/*  this function first computes the self-consistent temperature
 *  and abundance ratios, and then it calculates 
 *  (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRateFromU(double u, double rho, double *ne_guess)
{
    double temp;

    temp = convert_u_to_temp(u, rho, ne_guess);

    return CoolingRate(log10(temp), rho, ne_guess);
}

extern FILE *fd;


/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRate(double logT, double rho, double *nelec)
{
    double Lambda, Heat;
    double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
    double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
    double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
    double redshift;
    double T;

    if(logT <= Tmin)
        logT = Tmin + 0.5 * deltaT;	/* floor at Tmin */


    nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

    redshift = 1 / Time - 1;

    if(logT < Tmax)
    {
        find_abundances_and_rates(logT, rho, nelec);

        /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
        T = pow(10.0, logT);

        LambdaExcH0 = bH0 * ne * nH0;
        LambdaExcHep = bHep * ne * nHep;
        LambdaExc = LambdaExcH0 + LambdaExcHep;	/* excitation */

        LambdaIonH0 = 2.18e-11 * geH0 * ne * nH0;
        LambdaIonHe0 = 3.94e-11 * geHe0 * ne * nHe0;
        LambdaIonHep = 8.72e-11 * geHep * ne * nHep;
        LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* ionization */

        LambdaRecHp = 1.036e-16 * T * ne * (aHp * nHp);
        LambdaRecHep = 1.036e-16 * T * ne * (aHep * nHep);
        LambdaRecHepp = 1.036e-16 * T * ne * (aHepp * nHepp);
        LambdaRecHepd = 6.526e-11 * ad * ne * nHep;
        LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;

        LambdaFF = bff * (nHp + nHep + 4 * nHepp) * ne;

        Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

        LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / nHcgs;

        Lambda += LambdaCmptn;

        Heat = 0;
        if(J_UV != 0)
            Heat += (nH0 * epsH0 + nHe0 * epsHe0 + nHep * epsHep) / nHcgs;

    }
    else				/* here we're outside of tabulated rates, T>Tmax K */
    {
        /* at high T (fully ionized); only free-free and Compton cooling are present.  
           Assumes no heating. */

        Heat = 0;

        LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep =
            LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

        /* very hot: H and He both fully ionized */
        nHp = 1.0;
        nHep = 0;
        nHepp = yhelium;
        ne = nHp + 2.0 * nHepp;
        *nelec = ne;		/* note: in units of the hydrogen number density */

        T = pow(10.0, logT);
        LambdaFF =
            1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * ne;

        /* add inverse Compton cooling off the microwave background */
        LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / nHcgs;

        Lambda = LambdaFF + LambdaCmptn;
    }

    return (Heat - Lambda);
}





#define mymalloc(name, size) g_malloc(size)
void InitCoolMemory(void)
{
    BetaH0 = (double *) mymalloc("BetaH0", (NCOOLTAB + 1) * sizeof(double));
    BetaHep = (double *) mymalloc("BetaHep", (NCOOLTAB + 1) * sizeof(double));
    AlphaHp = (double *) mymalloc("AlphaHp", (NCOOLTAB + 1) * sizeof(double));
    AlphaHep = (double *) mymalloc("AlphaHep", (NCOOLTAB + 1) * sizeof(double));
    Alphad = (double *) mymalloc("Alphad", (NCOOLTAB + 1) * sizeof(double));
    AlphaHepp = (double *) mymalloc("AlphaHepp", (NCOOLTAB + 1) * sizeof(double));
    GammaeH0 = (double *) mymalloc("GammaeH0", (NCOOLTAB + 1) * sizeof(double));
    GammaeHe0 = (double *) mymalloc("GammaeHe0", (NCOOLTAB + 1) * sizeof(double));
    GammaeHep = (double *) mymalloc("GammaeHep", (NCOOLTAB + 1) * sizeof(double));
    Betaff = (double *) mymalloc("Betaff", (NCOOLTAB + 1) * sizeof(double));
}


void MakeCoolingTable(void)
    /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19 
       Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
    int i;
    double T;
    double Tfact;

#ifdef NEW_RATES
    double dE, P, A, X, K, U, T_eV;
    double b0, b1, b2, b3, b4, b5, c0, c1, c2, c3, c4, c5, y;	/* used in Scholz-Walter fit */
    double E1s_2, Gamma1s_2s, Gamma1s_2p;
#endif

    XH = 0.76;
    yhelium = (1 - XH) / (4 * XH);

    mhboltz = PROTONMASS / BOLTZMANN;

    if(All.MinGasTemp > 0.0)
        Tmin = log10(0.1 * All.MinGasTemp);
    else
        Tmin = 1.0;

    deltaT = (Tmax - Tmin) / NCOOLTAB;

    ethmin = pow(10.0, Tmin) * (1. + yhelium) / ((1. + 4. * yhelium) * mhboltz * GAMMA_MINUS1);
    /* minimum internal energy for neutral gas */

    for(i = 0; i <= NCOOLTAB; i++)
    {
        BetaH0[i] =
            BetaHep[i] =
            Betaff[i] =
            AlphaHp[i] = AlphaHep[i] = AlphaHepp[i] = Alphad[i] = GammaeH0[i] = GammaeHe0[i] = GammaeHep[i] = 0;


        T = pow(10.0, Tmin + deltaT * i);

        Tfact = 1.0 / (1 + sqrt(T / 1.0e5));

        if(118348 / T < 70)
            BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;

#ifdef NEW_RATES
        /* Scholtz-Walters 91 fit */
        if(T >= 2.0e3 && T < 1e8)
        {

            if(T >= 2.0e3 && T < 6.0e4)
            {
                b0 = -3.299613e1;
                b1 = 1.858848e1;
                b2 = -6.052265;
                b3 = 8.603783e-1;
                b4 = -5.717760e-2;
                b5 = 1.451330e-3;

                c0 = -1.630155e2;
                c1 = 8.795711e1;
                c2 = -2.057117e1;
                c3 = 2.359573;
                c4 = -1.339059e-1;
                c5 = 3.021507e-3;
            }
            else
            {
                if(T >= 6.0e4 && T < 6.0e6)
                {
                    b0 = 2.869759e2;
                    b1 = -1.077956e2;
                    b2 = 1.524107e1;
                    b3 = -1.080538;
                    b4 = 3.836975e-2;
                    b5 = -5.467273e-4;

                    c0 = 5.279996e2;
                    c1 = -1.939399e2;
                    c2 = 2.718982e1;
                    c3 = -1.883399;
                    c4 = 6.462462e-2;
                    c5 = -8.811076e-4;
                }
                else
                {
                    b0 = -2.7604708e3;
                    b1 = 7.9339351e2;
                    b2 = -9.1198462e1;
                    b3 = 5.1993362;
                    b4 = -1.4685343e-1;
                    b5 = 1.6404093e-3;

                    c0 = -2.8133632e3;
                    c1 = 8.1509685e2;
                    c2 = -9.4418414e1;
                    c3 = 5.4280565;
                    c4 = -1.5467120e-1;
                    c5 = 1.7439112e-3;
                }

                y = log(T);
                E1s_2 = 10.2;	/* eV */

                Gamma1s_2s =
                    exp(b0 + b1 * y + b2 * y * y + b3 * y * y * y + b4 * y * y * y * y + b5 * y * y * y * y * y);
                Gamma1s_2p =
                    exp(c0 + c1 * y + c2 * y * y + c3 * y * y * y + c4 * y * y * y * y + c5 * y * y * y * y * y);

                T_eV = T / eV_to_K;

                BetaH0[i] = E1s_2 * eV_to_erg * (Gamma1s_2s + Gamma1s_2p) * exp(-E1s_2 / T_eV);
            }
        }
#endif


        if(473638 / T < 70)
            BetaHep[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;

        Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));


#ifdef NEW_RATES
        AlphaHp[i] = 6.28e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);
#else
        AlphaHp[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);	/* old Cen92 fit */
#endif


        AlphaHep[i] = 1.5e-10 * pow(T, -0.6353);


#ifdef NEW_RATES
        AlphaHepp[i] = 3.36e-10 * pow(T / 1000, -0.2) / (1. + pow(T / 4.0e6, 0.7)) / sqrt(T);
#else
        AlphaHepp[i] = 4. * AlphaHp[i];	/* old Cen92 fit */
#endif

        if(470000 / T < 70)
            Alphad[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));


#ifdef NEW_RATES
        T_eV = T / eV_to_K;

        /* Voronov 97 fit */
        /* hydrogen */
        dE = 13.6;
        P = 0.0;
        A = 0.291e-7;
        X = 0.232;
        K = 0.39;

        U = dE / T_eV;
        GammaeH0[i] = A * (1.0 + P * sqrt(U)) * pow(U, K) * exp(-U) / (X + U);

        /* Helium */
        dE = 24.6;
        P = 0.0;
        A = 0.175e-7;
        X = 0.18;
        K = 0.35;

        U = dE / T_eV;
        GammaeHe0[i] = A * (1.0 + P * sqrt(U)) * pow(U, K) * exp(-U) / (X + U);

        /* Hellium II */
        dE = 54.4;
        P = 1.0;
        A = 0.205e-8;
        X = 0.265;
        K = 0.25;

        U = dE / T_eV;
        GammaeHep[i] = A * (1.0 + P * sqrt(U)) * pow(U, K) * exp(-U) / (X + U);

#else
        if(157809.1 / T < 70)
            GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;

        if(285335.4 / T < 70)
            GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;

        if(631515.0 / T < 70)
            GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
#endif

    }


}





/* table input (from file TREECOOL) for ionizing parameters */

#define JAMPL	1.0		/* amplitude factor relative to input table */
#define TABLESIZE 200		/* Max # of lines in TREECOOL */

static float inlogz[TABLESIZE];
static float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
static float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
static int nheattab;		/* length of table */


static double table[][7];

void ReadIonizeParams(char *fname)
{
    int i;
    for(i = 0; i < TABLESIZE; i++)
        gH0[i] = 0;


    if(fname == NULL) {
        for(i = 0; i < TABLESIZE; i++) {
            inlogz[i] = table[i][0];
            gH0[i] = table[i][1]; 
            gHe[i] = table[i][2];
            gHep[i] = table[i][3];
            eH0[i] = table[i][4];
            eHe[i] = table[i][5];
            eHep[i] = table[i][6];
            if(table[i][1] == 0.0) break;
        }
    } else {
        ROOTONLY {
            FILE *fdcool;

            if(!(fdcool = fopen(fname, "r")))
            {
                g_error(" Cannot read ionization table in file `%s'\n", fname);
            }

            for(i = 0; i < TABLESIZE; i++)
                if(fscanf(fdcool, "%g %g %g %g %g %g %g",
                            &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) != 7)
                    break;

            fclose(fdcool);
        }
        MPI_Bcast(inlogz, TABLESIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(gH0, TABLESIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(gHe, TABLESIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(gHep, TABLESIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(eHe, TABLESIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(eHep, TABLESIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    /*  nheattab is the number of entries in the table */

    for(i = 0, nheattab = 0; i < TABLESIZE; i++)
        if(gH0[i] != 0.0)
            nheattab++;
        else
            break;

    ROOTONLY {
        g_message("read ionization table with %d entries in file `%s'.", 
            nheattab, fname?fname:"<memory>");
    }
}

void LoadIonizeParams(double time)
{
    int i, ilow;
    double logz, dzlow, dzhi;
    double redshift;
    if (time == Time) return;
    Time = time;
    redshift = 1 / time - 1;

    logz = log10(redshift + 1.0);
    ilow = 0;
    for(i = 0; i < nheattab; i++)
    {
        if(inlogz[i] < logz)
            ilow = i;
        else
            break;
    }

    dzlow = logz - inlogz[ilow];
    dzhi = inlogz[ilow + 1] - logz;

    if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
        gJHe0 = gJHep = gJH0 = 0;
        epsHe0 = epsHep = epsH0 = 0;
        J_UV = 0;
        return;
    }
    else
        J_UV = 1.e-21;		/* irrelevant as long as it's not 0 */

    gJH0 = JAMPL * pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ilow + 1])) / (dzlow + dzhi));
    gJHe0 = JAMPL * pow(10., (dzhi * log10(gHe[ilow]) + dzlow * log10(gHe[ilow + 1])) / (dzlow + dzhi));
    gJHep = JAMPL * pow(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ilow + 1])) / (dzlow + dzhi));
    epsH0 = JAMPL * pow(10., (dzhi * log10(eH0[ilow]) + dzlow * log10(eH0[ilow + 1])) / (dzlow + dzhi));
    epsHe0 = JAMPL * pow(10., (dzhi * log10(eHe[ilow]) + dzlow * log10(eHe[ilow + 1])) / (dzlow + dzhi));
    epsHep = JAMPL * pow(10., (dzhi * log10(eHep[ilow]) + dzlow * log10(eHep[ilow + 1])) / (dzlow + dzhi));

    return;
}


void SetZeroIonization(void)
{
    gJHe0 = gJHep = gJH0 = 0;
    epsHe0 = epsHep = epsH0 = 0;
    J_UV = 0;
    Time = -1;
}


void InitCool(char * filename)
{

    InitCoolMemory();
    MakeCoolingTable();

    ReadIonizeParams(filename);

    SetZeroIonization();
}


static double table[][7] = {
  {0.000,3.03516e-14,1.37296e-14,3.04873e-16,1.74434e-25,1.76233e-25,1.00198e-26},
  {0.005,3.20557e-14,1.47386e-14,3.14717e-16,1.85463e-25,1.87090e-25,1.03701e-26},
  {0.010,3.37379e-14,1.57232e-14,3.24518e-16,1.96306e-25,1.97710e-25,1.07195e-26},
  {0.015,3.54076e-14,1.66914e-14,3.34310e-16,2.07032e-25,2.08173e-25,1.10691e-26},
  {0.020,3.70746e-14,1.76519e-14,3.44133e-16,2.17717e-25,2.18565e-25,1.14202e-26},
  {0.025,3.87497e-14,1.86137e-14,3.54027e-16,2.28440e-25,2.28979e-25,1.17740e-26},
  {0.030,4.04442e-14,1.95867e-14,3.64036e-16,2.39288e-25,2.39514e-25,1.21319e-26},
  {0.035,4.21704e-14,2.05814e-14,3.74207e-16,2.50352e-25,2.50277e-25,1.24953e-26},
  {0.040,4.39415e-14,2.16092e-14,3.84589e-16,2.61733e-25,2.61381e-25,1.28660e-26},
  {0.045,4.57713e-14,2.26820e-14,3.95237e-16,2.73534e-25,2.72948e-25,1.32455e-26},
  {0.050,4.76748e-14,2.38127e-14,4.06206e-16,2.85868e-25,2.85108e-25,1.36356e-26},
  {0.055,4.96649e-14,2.50126e-14,4.17547e-16,2.98833e-25,2.97975e-25,1.40379e-26},
  {0.060,5.17433e-14,2.62842e-14,4.29269e-16,3.12445e-25,3.11577e-25,1.44526e-26},
  {0.065,5.39090e-14,2.76279e-14,4.41369e-16,3.26695e-25,3.25918e-25,1.48795e-26},
  {0.070,5.61603e-14,2.90436e-14,4.53847e-16,3.41576e-25,3.41001e-25,1.53185e-26},
  {0.075,5.84958e-14,3.05314e-14,4.66697e-16,3.57078e-25,3.56831e-25,1.57694e-26},
  {0.080,6.09147e-14,3.20921e-14,4.79923e-16,3.73198e-25,3.73418e-25,1.62321e-26},
  {0.085,6.34210e-14,3.37301e-14,4.93541e-16,3.89966e-25,3.90806e-25,1.67070e-26},
  {0.090,6.60202e-14,3.54508e-14,5.07577e-16,4.07420e-25,4.09053e-25,1.71947e-26},
  {0.095,6.87178e-14,3.72597e-14,5.22055e-16,4.25603e-25,4.28216e-25,1.76960e-26},
  {0.100,7.15197e-14,3.91629e-14,5.37001e-16,4.44557e-25,4.48359e-25,1.82115e-26},
  {0.105,7.44316e-14,4.11660e-14,5.52443e-16,4.64323e-25,4.69541e-25,1.87419e-26},
  {0.110,7.74567e-14,4.32728e-14,5.68396e-16,4.84924e-25,4.91808e-25,1.92875e-26},
  {0.115,8.05977e-14,4.54866e-14,5.84877e-16,5.06376e-25,5.15205e-25,1.98486e-26},
  {0.120,8.38574e-14,4.78108e-14,6.01902e-16,5.28699e-25,5.39777e-25,2.04256e-26},
  {0.125,8.72386e-14,5.02490e-14,6.19488e-16,5.51911e-25,5.65571e-25,2.10187e-26},
  {0.130,9.07449e-14,5.28052e-14,6.37654e-16,5.76035e-25,5.92640e-25,2.16283e-26},
  {0.135,9.43816e-14,5.54857e-14,6.56428e-16,6.01112e-25,6.21054e-25,2.22550e-26},
  {0.140,9.81554e-14,5.82973e-14,6.75838e-16,6.27188e-25,6.50889e-25,2.28995e-26},
  {0.145,1.02073e-13,6.12473e-14,6.95918e-16,6.54310e-25,6.82228e-25,2.35627e-26},
  {0.150,1.06141e-13,6.43433e-14,7.16699e-16,6.82530e-25,7.15153e-25,2.42451e-26},
  {0.155,1.10366e-13,6.75924e-14,7.38214e-16,7.11896e-25,7.49749e-25,2.49477e-26},
  {0.160,1.14754e-13,7.09995e-14,7.60494e-16,7.42442e-25,7.86090e-25,2.56710e-26},
  {0.165,1.19310e-13,7.45687e-14,7.83569e-16,7.74196e-25,8.24246e-25,2.64158e-26},
  {0.170,1.24039e-13,7.83043e-14,8.07471e-16,8.07192e-25,8.64293e-25,2.71828e-26},
  {0.175,1.28947e-13,8.22108e-14,8.32231e-16,8.41460e-25,9.06308e-25,2.79726e-26},
  {0.180,1.34040e-13,8.62930e-14,8.57885e-16,8.77036e-25,9.50374e-25,2.87861e-26},
  {0.185,1.39325e-13,9.05574e-14,8.84475e-16,9.13966e-25,9.96583e-25,2.96241e-26},
  {0.190,1.44808e-13,9.50112e-14,9.12045e-16,9.52302e-25,1.04503e-24,3.04876e-26},
  {0.195,1.50497e-13,9.96618e-14,9.40641e-16,9.92097e-25,1.09583e-24,3.13776e-26},
  {0.200,1.56400e-13,1.04517e-13,9.70311e-16,1.03341e-24,1.14907e-24,3.22951e-26},
  {0.205,1.62524e-13,1.09584e-13,1.00111e-15,1.07629e-24,1.20488e-24,3.32411e-26},
  {0.210,1.68878e-13,1.14873e-13,1.03309e-15,1.12079e-24,1.26336e-24,3.42171e-26},
  {0.215,1.75471e-13,1.20392e-13,1.06634e-15,1.16696e-24,1.32463e-24,3.52245e-26},
  {0.220,1.82312e-13,1.26150e-13,1.10090e-15,1.21487e-24,1.38882e-24,3.62649e-26},
  {0.225,1.89411e-13,1.32158e-13,1.13687e-15,1.26456e-24,1.45604e-24,3.73398e-26},
  {0.230,1.96778e-13,1.38425e-13,1.17430e-15,1.31609e-24,1.52642e-24,3.84509e-26},
  {0.235,2.04420e-13,1.44957e-13,1.21328e-15,1.36952e-24,1.60010e-24,3.95998e-26},
  {0.240,2.12348e-13,1.51762e-13,1.25387e-15,1.42491e-24,1.67719e-24,4.07880e-26},
  {0.245,2.20569e-13,1.58846e-13,1.29613e-15,1.48232e-24,1.75781e-24,4.20173e-26},
  {0.250,2.29093e-13,1.66216e-13,1.34015e-15,1.54182e-24,1.84212e-24,4.32893e-26},
  {0.255,2.37931e-13,1.73879e-13,1.38603e-15,1.60347e-24,1.93024e-24,4.46063e-26},
  {0.260,2.47094e-13,1.81844e-13,1.43393e-15,1.66735e-24,2.02234e-24,4.59720e-26},
  {0.265,2.56597e-13,1.90119e-13,1.48409e-15,1.73355e-24,2.11859e-24,4.73907e-26},
  {0.270,2.66456e-13,1.98715e-13,1.53671e-15,1.80214e-24,2.21919e-24,4.88670e-26},
  {0.275,2.76684e-13,2.07640e-13,1.59204e-15,1.87323e-24,2.32432e-24,5.04058e-26},
  {0.280,2.87295e-13,2.16903e-13,1.65028e-15,1.94689e-24,2.43415e-24,5.20114e-26},
  {0.285,2.98297e-13,2.26504e-13,1.71157e-15,2.02315e-24,2.54877e-24,5.36861e-26},
  {0.290,3.09693e-13,2.36443e-13,1.77598e-15,2.10204e-24,2.66824e-24,5.54317e-26},
  {0.295,3.21489e-13,2.46719e-13,1.84363e-15,2.18357e-24,2.79263e-24,5.72499e-26},
  {0.300,3.33687e-13,2.57331e-13,1.91459e-15,2.26777e-24,2.92200e-24,5.91426e-26},
  {0.305,3.46292e-13,2.68275e-13,1.98896e-15,2.35463e-24,3.05640e-24,6.11109e-26},
  {0.310,3.59296e-13,2.79549e-13,2.06668e-15,2.44414e-24,3.19585e-24,6.31526e-26},
  {0.315,3.72693e-13,2.91149e-13,2.14768e-15,2.53628e-24,3.34037e-24,6.52651e-26},
  {0.320,3.86474e-13,3.03070e-13,2.23188e-15,2.63101e-24,3.48997e-24,6.74449e-26},
  {0.325,4.00628e-13,3.15309e-13,2.31918e-15,2.72828e-24,3.64466e-24,6.96887e-26},
  {0.330,4.15148e-13,3.27859e-13,2.40950e-15,2.82807e-24,3.80444e-24,7.19932e-26},
  {0.335,4.30037e-13,3.40716e-13,2.50279e-15,2.93035e-24,3.96932e-24,7.43572e-26},
  {0.340,4.45297e-13,3.53874e-13,2.59903e-15,3.03515e-24,4.13931e-24,7.67797e-26},
  {0.345,4.60935e-13,3.67327e-13,2.69817e-15,3.14243e-24,4.31441e-24,7.92600e-26},
  {0.350,4.76953e-13,3.81068e-13,2.80019e-15,3.25221e-24,4.49461e-24,8.17967e-26},
  {0.355,4.93350e-13,3.95087e-13,2.90502e-15,3.36442e-24,4.67986e-24,8.43878e-26},
  {0.360,5.10092e-13,4.09358e-13,3.01249e-15,3.47887e-24,4.86986e-24,8.70267e-26},
  {0.365,5.27139e-13,4.23848e-13,3.12243e-15,3.59530e-24,5.06427e-24,8.97056e-26},
  {0.370,5.44447e-13,4.38524e-13,3.23464e-15,3.71346e-24,5.26271e-24,9.24160e-26},
  {0.375,5.61966e-13,4.53349e-13,3.34890e-15,3.83305e-24,5.46475e-24,9.51489e-26},
  {0.380,5.79656e-13,4.68292e-13,3.46479e-15,3.95381e-24,5.67002e-24,9.78927e-26},
  {0.385,5.97501e-13,4.83345e-13,3.58126e-15,4.07565e-24,5.87839e-24,1.00628e-25},
  {0.390,6.15495e-13,4.98501e-13,3.69703e-15,4.19851e-24,6.08976e-24,1.03333e-25},
  {0.395,6.33628e-13,5.13759e-13,3.81073e-15,4.32233e-24,6.30403e-24,1.05985e-25},
  {0.400,6.51892e-13,5.29112e-13,3.92091e-15,4.44704e-24,6.52109e-24,1.08557e-25},
  {0.405,6.70266e-13,5.44548e-13,4.02635e-15,4.57250e-24,6.74074e-24,1.11029e-25},
  {0.410,6.88697e-13,5.60022e-13,4.12706e-15,4.69834e-24,6.96240e-24,1.13388e-25},
  {0.415,7.07118e-13,5.75479e-13,4.22336e-15,4.82413e-24,7.18542e-24,1.15628e-25},
  {0.420,7.25457e-13,5.90861e-13,4.31558e-15,4.94937e-24,7.40907e-24,1.17741e-25},
  {0.425,7.43639e-13,6.06104e-13,4.40409e-15,5.07356e-24,7.63257e-24,1.19718e-25},
  {0.430,7.61593e-13,6.21150e-13,4.48926e-15,5.19622e-24,7.85520e-24,1.21555e-25},
  {0.435,7.79283e-13,6.35966e-13,4.57121e-15,5.31709e-24,8.07643e-24,1.23264e-25},
  {0.440,7.96680e-13,6.50521e-13,4.65009e-15,5.43595e-24,8.29579e-24,1.24861e-25},
  {0.445,8.13752e-13,6.64782e-13,4.72602e-15,5.55253e-24,8.51278e-24,1.26364e-25},
  {0.450,8.30465e-13,6.78717e-13,4.79916e-15,5.66661e-24,8.72686e-24,1.27789e-25},
  {0.455,8.46779e-13,6.92287e-13,4.86956e-15,5.77787e-24,8.93742e-24,1.29153e-25},
  {0.460,8.62630e-13,7.05442e-13,4.93682e-15,5.88589e-24,9.14370e-24,1.30445e-25},
  {0.465,8.77947e-13,7.18124e-13,5.00044e-15,5.99020e-24,9.34488e-24,1.31651e-25},
  {0.470,8.92654e-13,7.30275e-13,5.05989e-15,6.09029e-24,9.54006e-24,1.32757e-25},
  {0.475,9.06670e-13,7.41833e-13,5.11460e-15,6.18564e-24,9.72833e-24,1.33747e-25},
  {0.480,9.19911e-13,7.52732e-13,5.16403e-15,6.27569e-24,9.90870e-24,1.34606e-25},
  {0.485,9.32296e-13,7.62912e-13,5.20780e-15,6.35990e-24,1.00802e-23,1.35328e-25},
  {0.490,9.43739e-13,7.72311e-13,5.24556e-15,6.43770e-24,1.02419e-23,1.35906e-25},
  {0.495,9.54150e-13,7.80863e-13,5.27693e-15,6.50848e-24,1.03926e-23,1.36336e-25},
  {0.500,9.63436e-13,7.88500e-13,5.30153e-15,6.57161e-24,1.05313e-23,1.36613e-25},
  {0.505,9.71506e-13,7.95155e-13,5.31900e-15,6.62650e-24,1.06569e-23,1.36730e-25},
  {0.510,9.78299e-13,8.00783e-13,5.32906e-15,6.67274e-24,1.07686e-23,1.36681e-25},
  {0.515,9.83755e-13,8.05339e-13,5.33145e-15,6.70994e-24,1.08658e-23,1.36463e-25},
  {0.520,9.87817e-13,8.08778e-13,5.32590e-15,6.73772e-24,1.09477e-23,1.36068e-25},
  {0.525,9.90421e-13,8.11053e-13,5.31214e-15,6.75566e-24,1.10135e-23,1.35491e-25},
  {0.530,9.91509e-13,8.12119e-13,5.28993e-15,6.76339e-24,1.10625e-23,1.34728e-25},
  {0.535,9.91034e-13,8.11944e-13,5.25917e-15,6.76061e-24,1.10942e-23,1.33775e-25},
  {0.540,9.88954e-13,8.10498e-13,5.21981e-15,6.74704e-24,1.11081e-23,1.32634e-25},
  {0.545,9.85225e-13,8.07751e-13,5.17181e-15,6.72241e-24,1.11036e-23,1.31305e-25},
  {0.550,9.79804e-13,8.03673e-13,5.11512e-15,6.68643e-24,1.10803e-23,1.29786e-25},
  {0.555,9.72656e-13,7.98241e-13,5.04976e-15,6.63889e-24,1.10377e-23,1.28080e-25},
  {0.560,9.63786e-13,7.91453e-13,4.97591e-15,6.57980e-24,1.09756e-23,1.26189e-25},
  {0.565,9.53210e-13,7.83316e-13,4.89379e-15,6.50925e-24,1.08941e-23,1.24117e-25},
  {0.570,9.40948e-13,7.73839e-13,4.80367e-15,6.42734e-24,1.07932e-23,1.21870e-25},
  {0.575,9.27023e-13,7.63033e-13,4.70584e-15,6.33420e-24,1.06729e-23,1.19453e-25},
  {0.580,9.11479e-13,7.50925e-13,4.60072e-15,6.23011e-24,1.05334e-23,1.16874e-25},
  {0.585,8.94420e-13,7.37592e-13,4.48914e-15,6.11574e-24,1.03758e-23,1.14153e-25},
  {0.590,8.75971e-13,7.23127e-13,4.37207e-15,5.99193e-24,1.02012e-23,1.11313e-25},
  {0.595,8.56270e-13,7.07630e-13,4.25056e-15,5.85960e-24,1.00109e-23,1.08379e-25},
  {0.600,8.35464e-13,6.91213e-13,4.12575e-15,5.71971e-24,9.80648e-24,1.05377e-25},
  {0.605,8.13689e-13,6.73973e-13,3.99881e-15,5.57316e-24,9.58915e-24,1.02335e-25},
  {0.610,7.91014e-13,6.55949e-13,3.87085e-15,5.42035e-24,9.35945e-24,9.92793e-26},
  {0.615,7.67495e-13,6.37169e-13,3.74305e-15,5.26160e-24,9.11776e-24,9.62355e-26},
  {0.620,7.43195e-13,6.17664e-13,3.61667e-15,5.09725e-24,8.86454e-24,9.32324e-26},
  {0.625,7.18183e-13,5.97467e-13,3.49305e-15,4.92770e-24,8.60025e-24,9.03005e-26},
  {0.630,6.92510e-13,5.76609e-13,3.37306e-15,4.75324e-24,8.32529e-24,8.74586e-26},
  {0.635,6.66160e-13,5.55097e-13,3.25568e-15,4.57383e-24,8.03966e-24,8.46801e-26},
  {0.640,6.39103e-13,5.32935e-13,3.13941e-15,4.38935e-24,7.74325e-24,8.19269e-26},
  {0.645,6.11305e-13,5.10128e-13,3.02264e-15,4.19968e-24,7.43598e-24,7.91577e-26},
  {0.650,5.82733e-13,4.86680e-13,2.90365e-15,4.00470e-24,7.11776e-24,7.63286e-26},
  {0.655,5.53398e-13,4.62628e-13,2.78116e-15,3.80458e-24,6.78897e-24,7.34068e-26},
  {0.660,5.23475e-13,4.38110e-13,2.65588e-15,3.60049e-24,6.45155e-24,7.04071e-26},
  {0.665,4.93190e-13,4.13298e-13,2.52895e-15,3.39394e-24,6.10796e-24,6.73562e-26},
  {0.670,4.62781e-13,3.88374e-13,2.40165e-15,3.18652e-24,5.76080e-24,6.42827e-26},
  {0.675,4.32506e-13,3.63534e-13,2.27531e-15,2.97993e-24,5.41289e-24,6.12171e-26},
  {0.680,4.02608e-13,3.38965e-13,2.15110e-15,2.77581e-24,5.06689e-24,5.81862e-26},
  {0.685,3.73229e-13,3.14778e-13,2.02940e-15,2.57509e-24,4.72449e-24,5.51973e-26},
  {0.690,3.44498e-13,2.91080e-13,1.91041e-15,2.37866e-24,4.38723e-24,5.22542e-26},
  {0.695,3.16553e-13,2.67978e-13,1.79434e-15,2.18744e-24,4.05678e-24,4.93605e-26},
  {0.700,2.89539e-13,2.45592e-13,1.68142e-15,2.00239e-24,3.73490e-24,4.65201e-26},
  {0.705,2.63589e-13,2.24030e-13,1.57187e-15,1.82444e-24,3.42325e-24,4.37371e-26},
  {0.710,2.38785e-13,2.03361e-13,1.46587e-15,1.65414e-24,3.12293e-24,4.10162e-26},
  {0.715,2.15194e-13,1.83648e-13,1.36361e-15,1.49198e-24,2.83496e-24,3.83623e-26},
  {0.720,1.92888e-13,1.64957e-13,1.26529e-15,1.33849e-24,2.56041e-24,3.57804e-26},
  {0.725,1.71944e-13,1.47356e-13,1.17110e-15,1.19420e-24,2.30039e-24,3.32758e-26},
  {0.730,1.52421e-13,1.30901e-13,1.08119e-15,1.05956e-24,2.05587e-24,3.08532e-26},
  {0.735,1.34328e-13,1.15602e-13,9.95540e-16,9.34625e-25,1.82712e-24,2.85142e-26},
  {0.740,1.17654e-13,1.01458e-13,9.14061e-16,8.19349e-25,1.61426e-24,2.62600e-26},
  {0.745,1.02391e-13,8.84631e-14,8.36666e-16,7.13671e-25,1.41738e-24,2.40916e-26},
  {0.750,8.85234e-14,7.66120e-14,7.63247e-16,6.17513e-25,1.23658e-24,2.20098e-26},
  {0.755,7.60273e-14,6.58897e-14,6.93737e-16,5.30720e-25,1.07181e-24,2.00163e-26},
  {0.760,6.48482e-14,5.62575e-14,6.28219e-16,4.52939e-25,9.22646e-25,1.81158e-26},
  {0.765,5.49201e-14,4.76680e-14,5.66814e-16,3.83742e-25,7.88556e-25,1.63134e-26},
  {0.770,4.61705e-14,4.00688e-14,5.09642e-16,3.22657e-25,6.68935e-25,1.46149e-26},
  {0.775,3.85199e-14,3.34019e-14,4.56825e-16,2.69166e-25,5.63114e-25,1.30260e-26},
  {0.780,3.18847e-14,2.76053e-14,4.08322e-16,2.22719e-25,4.70327e-25,1.15496e-26},
  {0.785,2.61848e-14,2.26161e-14,3.63550e-16,1.82784e-25,3.89653e-25,1.01791e-26},
  {0.790,2.13351e-14,1.83663e-14,3.21775e-16,1.48789e-25,3.20070e-25,8.90547e-27},
  {0.795,1.72430e-14,1.47819e-14,2.82204e-16,1.20115e-25,2.60465e-25,7.71863e-27},
  {0.800,1.38082e-14,1.17826e-14,2.43984e-16,9.60883e-26,2.09630e-25,6.60776e-27},
  {0.805,1.09348e-14,9.28974e-15,2.06562e-16,7.60575e-26,1.66418e-25,5.56549e-27},
  {0.810,8.55931e-15,7.24513e-15,1.70516e-16,5.95623e-26,1.30119e-25,4.59785e-27},
  {0.815,6.62143e-15,5.59106e-15,1.36700e-16,4.61541e-26,1.00067e-25,3.71403e-27},
  {0.820,5.05572e-15,4.26510e-15,1.06025e-16,3.53472e-26,7.55377e-26,2.92375e-27},
  {0.825,3.79149e-15,3.19987e-15,7.94636e-17,2.66172e-26,5.57487e-26,2.23726e-27},
  {0.830,2.76009e-15,2.33100e-15,5.76164e-17,1.94610e-26,3.98834e-26,1.65825e-27},
  {0.835,1.91144e-15,1.61533e-15,3.97497e-17,1.35327e-26,2.71614e-26,1.16780e-27},
  {0.840,1.19633e-15,1.01154e-15,2.47998e-17,8.49865e-27,1.67609e-26,7.41001e-28},
  {0.845,5.61822e-16,4.75209e-16,1.16223e-17,4.00000e-27,7.79785e-27,3.51094e-28},
  {0.850,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00},
  {0.0,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00},
};
