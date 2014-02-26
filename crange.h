/**
 * \file crange.h
 * \brief Header file for crange.
 *
 * This header file collects all the other header files needed to compile
 * crange, as well as all defines, function declarations, etc.
 */
/*
 * Standard includes here to keep crange.c clean
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h> /* to get getopt() */
#include <errno.h>
#include <time.h>
/*
 * Include headers for various complex arithmetic functions.
 */
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
/*
 * Include config.h
 */
#include <config.h>
/*
 * Include header for parsing ini files.
 */
#ifdef HAVE_INIPARSER_H
#include <iniparser.h>
#endif
/*
 * These values define the range-energy tables
 */
#define LOGTENEMIN 0.0 /**< \f$ \log_{10} E_{\mathrm{min}} \f$ Minimum energy in units of A MeV. */
#define LOGTENEMAX 6.0 /**< \f$ \log_{10} E_{\mathrm{max}} \f$ Maximum energy in units of A MeV. */
#define MAXE 200 /**< The number of energies in the range-energy tables. */
#define MAXAB 50 /**< The number of range tables.  Arbitrary, but should be larger than the number of targets. */
/*
 * M_PI, M_PI_2 and M_LN10 should be defined in math.h.
 */
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /**< The value of pi, in case it is not defined in math.h. */
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169163975144 /**< The value of pi/2, in case it is not defined in math.h. */
#endif
#ifndef M_LN10
#define M_LN10 2.30258509299404568402 /**< The value of \f$ \ln 10 \f$ in case it is not defined in math.h */
#endif
/*
 * Physical constants.
 */
#ifdef ALPHA
#undef ALPHA
#endif
#ifndef ALPHA
#define ALPHA 7.29735301383e-3 /**< The fine structure constant. */
#endif
#ifdef ATOMICMASSUNIT
#undef ATOMICMASSUNIT
#endif
#ifndef ATOMICMASSUNIT
#define ATOMICMASSUNIT 931.4943 /**< 1 amu in units of MeV/c<sup>2</sup>. */
#endif
#ifdef PROTONMASS
#undef PROTONMASS
#endif
#ifndef PROTONMASS
#define PROTONMASS 938.2723 /**< Proton mass in units of MeV/c<sup>2</sup>. */
#endif
#ifdef ELECTRONMASS
#undef ELECTRONMASS
#endif
#ifndef ELECTRONMASS
#define ELECTRONMASS 0.511003e+6 /**< Electron mass in units of eV/c<sup>2</sup>. */
#endif
/*
 * These define the values in the switch bit field.
 */
#define SSWITCH_BA  0x001 /**< Barkas effect bit. */
#define SSWITCH_SH  0x002 /**< Shell effect bit. */
#define SSWITCH_LE  0x004 /**< Leung effect bit. */
#define SSWITCH_ND  0x008 /**< Density effect switcher bit. */
#define SSWITCH_EC  0x010 /**< Electron capture switcher bit. */
#define SSWITCH_NS  0x020 /**< Finite Nuclear Size effect bit. */
#define SSWITCH_KI  0x040 /**< Kinetic effect bit. */
#define SSWITCH_RA  0x080 /**< Radiative correction effect bit. */
#define SSWITCH_PA  0x100 /**< Pair Production energy loss bit. */
#define SSWITCH_BR  0x200 /**< Projectile Bremsstrahlung bit. */
#define SSWITCH_DEFAULT (SSWITCH_ND | SSWITCH_NS) /**< Default bits set density effect and finite nuclear size. */
/*
 *
 */
#define NAMEWIDTH 8 /**< The maximum number of characters in a target name. */
/*
 *
 */
/**
 * \brief Structure containing target data.
 *
 * This structure contains all the data related to target materials.
 */
struct TDATA {
    /**
     * \name Material name
     * @{
     */
    char name[NAMEWIDTH+1]; /**< The name of the material. */
    /** @} */
    /**
     * \name General parameters
     * @{
     */
    double z2;   /**< The mean nuclear charge. */
    double a2;   /**< The mean atomic mass number. */
    double iadj; /**< The logarithmic mean ionization potential [eV]. */
    double rho;  /**< The density [g cm<sup>-3</sup>]. */
    double pla;  /**< The plasma frequency [eV]. */
    double etad; /**< The ratio of density to density at STP for gaseous targets. Should be set to zero for non-gaseous materials. */
    double bind; /**< The total electronic binding energy [keV]. */
    /** @} */
    /**
     * \name Density effect parameters
     * @{
     */
    double X0; /**< Value of \f$ \log_{10} \beta\gamma \f$ at which the density effect turns on. */
    double X1; /**< Value of \f$ \log_{10} \beta\gamma \f$ above which the high-energy form of the density effect may be used. */
    double a;  /**< Parameter used to interpolate the density effect between the values of X0 and X1. */
    double m;  /**< Parameter used to interpolate the density effect between the values of X0 and X1. */
    double d0; /**< Low-energy density effect parameter, only non-zero for conducting materials. */
    /** @} */
};
/**
 * \brief Define tdata.
 *
 * Define a tdata variable for convenience.
 *
 */
typedef struct TDATA tdata;
/**
 * \brief Structure to store range tables.
 *
 * This structure contains a range table and its associated metadata.
 */
struct RANGE_TABLE {
    double z1;          /**< The projectile charge. */
    double a1;          /**< The projectile mass. */
    short sswitch;      /**< The switch bit field. */
    tdata *target;      /**< A pointer to the ::TDATA structure used in constructing the table. */
    time_t timestamp;   /**< The time at which the table was created. */
    double range[MAXE]; /**< The actual table of range values. */
};
/**
 * \brief Define range_table.
 *
 * Define a range_table variable for convenience.
 */
typedef struct RANGE_TABLE range_table;
/**
 * \brief The range-energy table.
 *
 * This external variable contains the range-energy tables.
 */
range_table trange[MAXAB];
/*
 * Delcare all the functions in crange.c
 */
gsl_complex complex_hyperg( gsl_complex a, gsl_complex b, gsl_complex z );
gsl_complex complex_lngamma( gsl_complex z );
double effective_charge( double z0, double e1, double z2, short sswitch );
double djdx( double e1, double z0, double I0, double f0, double K, short sswitch, tdata *target);
double dedx( double e1, double rel0, double z0, double a1, short sswitch, tdata *target );
double delta( double g, tdata *target );
double olddelta( double g, tdata *target );
double bma( double z1, double b );
double relbloch( double z12, double b1, double lambda, double theta0 );
double lindhard( double zz, double aa, double bb, short sswitch );
double Fbrems( double x );
double range( double e, double z1, double a1, short sswitch, tdata *target, int *tno );
double qrange( double e, double z1, double a1, short sswitch, tdata *target );
double benton( double e, double z1, double a1, tdata *target );
double renergy( double e, double r0, double z1, double a1, short sswitch, tdata *target );
void run_range( FILE *finput, FILE *foutput, short sswitch, tdata *extratargets );
short init_switch( char *switchfile);
tdata *init_target( char *targetfile );
void init_table(void);
double energy_table( int i );
tdata *find_target( char *target, tdata *extratargets );
void print_target( tdata *target );
