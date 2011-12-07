/*
 * This is version 1.5.3 of the Berkeley Range-Energy Calculator
 *
 * Benjamin Weaver, weaver@SSL.Berkeley.EDU
 * Space Sciences Laboratory
 * University of California
 * Berkeley, CA 94720-7450
 * http://ultraman.ssl.berkeley.edu/~weaver/dedx/
 *
 * Copyright (C) 2001-2011 Benjamin Weaver
 *
 *
 * This program is free software which I release under the GNU Lesser
 * General Public License. You may redistribute and/or modify this program
 * under the terms of that license as published by the Free Software
 * Foundation; either version 2.1 of the License, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * To get a copy of the GNU Lesser General Puplic License, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
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
 * Define working directory.  This is the directory in which the
 * binary lives as well as the data files.
 * IMPORTANT:  Change this to your own working directory!
 */
#ifdef CRANGE_DIR
#undef CRANGE_DIR
#endif
#ifndef CRANGE_DIR
#define CRANGE_DIR "."
#endif
/*
 * MAXE sets the number of energies in the range-energy tables.
 */
#define MAXE 200
/*
 * MAXAB sets the number of different target media which can be stored
 * in the range-energy tables.  It is suggested that MAXAB > number of
 * target media in the target.dat file.
 */
#define MAXAB 50
/*
 * Miscellaneous defines. ALPHA is the fine structure constant, and we want
 * exactly that value.  Otherwise, M_PI and M_PI_2 should be defined in math.h.
 */
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169163975144
#endif
#ifdef ALPHA
#undef ALPHA
#endif
#ifndef ALPHA
#define ALPHA 7.29735301383e-3
#endif
/*
 * Define switch bits.
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
 * These external variables contain the range-energy tables.
 */
double trange[MAXE][MAXAB];
double tenerg[MAXE];
/**
 * @struct TDATA crange.h
 * @brief Structure containing target data.
 *
 * Structure containing target data.
 */
typedef struct TDATA {
    /**
     * @name General parameters.
     */
    /* @{ */
    double z2;   /**< The mean nuclear charge. */
    double a2;   /**< The mean atomic mass number. */
    double iadj; /**< The logarithmic mean ionization potential [eV]. */
    double rho;  /**< The density [g cm^-3]. */
    double pla;  /**< The plasma frequency [eV]. */
    double etad; /**< The ratio of density to density at STP for gaseous targets. Should be set to zero for non-gaseous materials. */
    double bind; /**< The total electronic binding energy [eV]. */
    /* @} */
    /**
     * @name Density effect parameters.
     */
    /* @{ */
    double X0; /**< Value of @f$\log_{10} \beta\gamma@f$ at which the density effect turns on. */
    double X1; /**< Value of @f$\log_{10} \beta\gamma@f$ above which the high-energy form of the density effect may be used. */
    double a;  /**< Parameter used to interpolate the density effect between the values of X0 and X1. */
    double m;  /**< Parameter used to interpolate the density effect between the values of X0 and X1. */
    double d0; /**< Low-energy density effect parameter, only non-zero for conducting materials. */
    /* @} */
    /**
     * @name Material name
     */
    /* @{ */
    char tname[9]; /**< The name of the material. */
    /* @} */
} tdata;
/*
 * Declare an external array of structures
 */
tdata t[MAXAB];
/*
 * Delcare all the functions in crange.c
 */
gsl_complex complex_hyperg( gsl_complex a, gsl_complex b, gsl_complex z );
gsl_complex complex_lngamma( gsl_complex z );
double effective_charge( double z0, double e1, double z2, short sswitch );
double djdx( double e1, double z0, double I0, double f0, double K, short sswitch, int tno);
double dedx( double e1, double rel0, double z0, double a1, short sswitch, int tno );
double delta( double g, int tno );
double olddelta( double g, int tno );
double bma( double z1, double b );
double relbloch( double z12, double b1, double lambda, double theta0 );
double lindhard( double zz, double aa, double bb, short sswitch );
double Fbrems( double x );
double range( double e, double z1, double a1, short sswitch, int tno );
double qrange( double e, double z1, double a1, short sswitch, int tno );
double benton( double e, double z1, double a1, int tno );
double renergy( double e, double r0, double z1, double a1, short sswitch, int tno );
void run_range( FILE *finput, FILE *foutput, short sswitch );
short init_switch( char *switchfile);
int init_tables( char *targetfile );
