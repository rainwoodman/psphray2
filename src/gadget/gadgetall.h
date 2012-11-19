
/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{

  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */

  /* system of units  */

  double UnitTime_in_s,		/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/*!< factor to convert internal mass unit to grams/h */
    UnitLength_in_cm,		/*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/*!< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitEnergy_in_cgs;		/*!< factor to convert internal energy to cgs units */

} All;

/* ... often used physical constants (cgs units) */
#define  HYDROGEN_MASSFRAC 0.76	/*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define MAXITER 150
#define  GAMMA         (5.0/3.0)	/*!< adiabatic index of simulated gas */
#define  GAMMA_MINUS1  (GAMMA-1)
static const double  GRAVITY    = 6.672e-8;
static const double  SOLAR_MASS = 1.989e33;
static const double  SOLAR_LUM  = 3.826e33;
static const double  RAD_CONST  = 7.565e-15;
static const double  AVOGADRO   = 6.0222e23;
static const double  BOLTZMANN  = 1.38066e-16;
static const double  GAS_CONST  = 8.31425e7;
static const double  C          = 2.9979e10;
static const double  PLANCK     = 6.6262e-27;
static const double  CM_PER_MPC = 3.085678e24;
static const double  PROTONMASS = 1.6726e-24;
static const double  ELECTRONMASS= 9.10953e-28;
static const double  THOMPSON    = 6.65245e-25;
static const double  ELECTRONCHARGE = 4.8032e-10;
static const double  HUBBLE         = 3.2407789e-18	/* in h/sec */;
static const double  LYMAN_ALPHA     = 1215.6e-8	/* 1215.6 Angstroem */;
static const double  LYMAN_ALPHA_HeII = 303.8e-8	/* 303.8 Angstroem */;
static const double  OSCILLATOR_STRENGTH      = 0.41615;
static const double  OSCILLATOR_STRENGTH_HeII = 0.41615;

double convert_u_to_temp(double u, double rho, double *ne_guess);
double CoolingRate(double logT, double rho, double *nelec);
double CoolingRateFromU(double u, double rho, double *ne_guess);
double GetCoolingTime(double u_old, double rho,  double *ne_guess);

void   find_abundances_and_rates(double logT, double rho, double *ne_guess);
void   InitCool(char * treecool);
void   LoadIonizeParams(double time);
void   SetZeroIonization(void);

