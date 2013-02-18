#include <glib.h>
#include <stdint.h>
#include <gsl/gsl_integration.h>
#include "commonblock.h"
#include <math.h>

static double R8;
static double r_tophat;

static double AA, BB, CC;
static double nu;
static double Norm;

static int NPowerTable;

static struct pow_table
{
  double logk, logD;
}
 *PowerTable;


static double PowerSpec_Tabulated(double k);
static double PowerSpec_EH(double k);
static double PowerSpec_Efstathiou(double k);
double GrowthFactor(double astart, double aend);
static double growth(double a);
static double growth_int(double a);
static double TopHatSigma2(double R);

double PowerSpec(double k)
{
  double power, alpha, Tf;
  switch (CB.IC.WhichSpectrum)
    {
    case 1:
      power = PowerSpec_EH(k);
      break;

    case 2:
      g_error("tabulated powerspec is unsupported");
      power = PowerSpec_Tabulated(k);
      break;

    default:
      power = PowerSpec_Efstathiou(k);
      break;
    }


  power *= pow(k, CB.IC.PrimordialIndex - 1.0);

  return power;
}

double GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}

int compare_logk(const void *a, const void *b)
{
  if(((struct pow_table *) a)->logk < (((struct pow_table *) b)->logk))
    return -1;

  if(((struct pow_table *) a)->logk > (((struct pow_table *) b)->logk))
    return +1;

  return 0;
}

void init_power(void)
{
  AA = 6.4 / CB.IC.ShapeGamma * (1000 * CB.U.KPC_h);
  BB = 3.0 / CB.IC.ShapeGamma * (1000 * CB.U.KPC_h);
  CC = 1.7 / CB.IC.ShapeGamma * (1000 * CB.U.KPC_h);
  nu = 1.13;

  R8 = 8 * (1000 * CB.U.KPC_h); /* 8 Mpc/h */

  Norm = 1.0;
  Norm = CB.C.Sigma8 * CB.C.Sigma8 / TopHatSigma2(R8);
}


#if 0
static double tk_eh(double k)		/* from Martin White */
{
  static int inited = 0;
  static double A, S;
  static double MPC_h;
  const double hubble = CB.C.h;
  const double omegam = CB.C.OmegaM;
  const double theta = 2.728 / 2.7;
  if(!inited) {
      /* other input parameters */

      const double omegab = CB.C.OmegaB;
      const double ommh2 = omegam * hubble * hubble;
      const double ombh2 = omegab * hubble * hubble;
      MPC_h = 1000 * CB.U.KPC_h;
      S = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
      A = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
        + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
      inited = 1;
  }

  /* convert to h/Mpc */
  const double ks = 0.43 * (k * MPC_h) * S;

  const double gamma = (A + (1. - A) / (1. + ks * ks * ks * ks)) * (omegam * hubble);
  const double q = (k * MPC_h) * theta * theta / gamma;
  const double L0 = log(2. * M_E + 1.8 * q);
  const double C0 = 14.2 + 731. / (1. + 62.5 * q);
  return L0 / (L0 + C0 * q * q);
}
#else
static double tk_eh(double k)		/* from Martin White */
{
  double q, theta, ommh2, a, s, gamma, L0, C0;
  double tmp;
  double omegam, ombh2, hubble;

  /* other input parameters */
  hubble = CB.C.h;

  omegam = CB.C.OmegaM;
  ombh2 = CB.C.OmegaB * hubble * hubble;

  k *= (1000 * CB.U.KPC_h);

  theta = 2.728 / 2.7;
  ommh2 = omegam * hubble * hubble;
  s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
  a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
    + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
  gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma *= omegam * hubble;
  q = k * theta * theta / gamma;
  L0 = log(2. * exp(1.) + 1.8 * q);
  C0 = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp);
}


#endif
static double PowerSpec_Tabulated(double k)
{
  double logk, logD, P, kold, u, dlogk, Delta2;
  int binlow, binhigh, binmid;

  double mydlogk,dlogk_PowerTable;
  int mybinhigh,mybinlow,mybinmid;

  kold = k;

#if 0
  This need to be fixed
  k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);	/* convert to h/Mpc */

#endif
  logk = log10(k);

  if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
    return 0;

  dlogk_PowerTable = PowerTable[1].logk-PowerTable[0].logk;
  mydlogk = logk - PowerTable[0].logk;
  mybinlow = (int)(mydlogk/dlogk_PowerTable);
  mybinhigh = mybinlow+1;

  dlogk = PowerTable[mybinhigh].logk - PowerTable[mybinlow].logk;

  if(dlogk == 0)
    g_error("dlogk == 0");

  u = (logk - PowerTable[mybinlow].logk) / dlogk;

  logD = (1 - u) * PowerTable[mybinlow].logD + u * PowerTable[mybinhigh].logD;

  P = Norm*pow(10.0, logD);//*2*M_PI*M_PI;

  return P;
}

static double PowerSpec_Efstathiou(double k)
{
  return Norm * k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}



static double PowerSpec_EH(double k)	/* Eisenstein & Hu */
{
  return Norm * k * pow(tk_eh(k), 2);
    
}



double sigma2_int(double k)
{
  double kr, kr3, kr2, w, x;

  kr = r_tophat * k;
  kr2 = kr * kr;
  kr3 = kr2 * kr;

  if(kr < 1e-8)
    return 0;

  w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
  x = 4 * G_PI * k * k * w * w * PowerSpec(k);

  return x;

}
static double qrombwrapper(double x, double (*func)(double x)) {
    return func(x);
}
static double qromb(double (*func)(double x), double a, double b) {
    gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = qrombwrapper;
    F.params = (void*) func;
    double result;
    double error;
    gsl_integration_qag(&F, a, b, 0, 1e-7, 1000, GSL_INTEG_GAUSS61, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);
    return result;
}

static double TopHatSigma2(double R)
{
  r_tophat = R;

  return qromb(sigma2_int, 0, 500.0 * 1 / R);	/* note: 500/R is here chosen as 
						   integration boundary (infinity) */
}



static double growth(double a)
{
  double hubble_a;

  hubble_a = sqrt(CB.C.OmegaM / (a * a * a) + (1 - CB.C.OmegaM - CB.C.OmegaL) / (a * a) + CB.C.OmegaL);

  return hubble_a * qromb(growth_int, 0, a);
}


static double growth_int(double a)
{
  return pow(a / (CB.C.OmegaM + (1 - CB.C.OmegaM - CB.C.OmegaL) * a + CB.C.OmegaL * a * a * a), 1.5);
}


double F_Omega(double a)
{
  double omega_a;

  omega_a = CB.C.OmegaM / (CB.C.OmegaM + a * (1 - CB.C.OmegaM - CB.C.OmegaL) + a * a * a * CB.C.OmegaL);

  return pow(omega_a, 0.6);
}

