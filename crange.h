/*
 * This is version 1.5.3 of the Berkeley Range-Energy Calculator
 *
 * Benjamin Weaver, weaver@SSL.Berkeley.EDU
 * Space Sciences Laboratory
 * University of California
 * Berkeley, CA 94720-7450
 * http://ultraman.ssl.berkeley.edu/~weaver/dedx/
 *
 * Copyright (C) 2001-2003 Benjamin Weaver
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
/*
 * Include headers for various complex arithmetic functions.
 */
#include <fcomplex.h>
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
#ifdef DEDX
#undef DEDX
#endif
#ifndef DEDX
#define DEDX "."
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
 * Miscellaneous defines. The ifdefs should insure that we get the
 * values listed here. ALPHA is the fine structure constant.
 */
#ifdef ALPHA
#undef ALPHA
#endif
#ifndef ALPHA
#define ALPHA 7.29735301383e-3
#endif

#define FBA  0x001 /* Barkas effect */
#define FSH  0x002 /* Shell effect */
#define FLE  0x004 /* Leung effect */
#define FD   0x008 /* Density effect switcher */
#define FE   0x010 /* Electron capture switcher */
#define FNS  0x020 /* Finite Nuclear Size effect */
#define FKIN 0x040 /* Kinetic effect */
#define FRAD 0x080 /* Radiative correction effect */
#define FPA  0x100 /* Pair Production energy loss */
#define FBR  0x200 /* Projectile Bremsstrahlung */
/*
 * External variable for holding calculation switches.
 */
short sswitch;
/*
 * These external variables contain the range-energy tables.
 */
double trange[MAXE][MAXAB];
double tenerg[MAXE];
/*
 * This structure holds the data on the targest read from the target.dat
 * file.
 * z2:     target material mean nuclear charge
 * a2:     target material mean atomic number
 * iadj:   target material logarithmic mean ionization potential (eV)
 * rho:    target material density in g cm^-3
 * pla:    target material plasma frequency in eV
 * etad:   target material ratio of density to density at STP (for gasses)
 *           target materials which are not gasses should have etad = 0
 * bind:   target material total electronic binding energy in eV
 * X0,X1,a,m,d0:  parameters for the density effect calculation
 * tname:   names of target materials
 */
typedef struct TDATA {
    double z2,a2,iadj,rho,pla,X0,X1,a,m,d0,etad,bind;
    char tname[9];
} tdata;
/*
 * Declare an external array of structures
 */
tdata t[MAXAB];
/*
 * Delcare all the functions in crange.c
 */
double dedx( double e1, double rel0, double z0, double a1, int tno );
double delta( double g, int tno );
double olddelta( double g, int tno );
double bma( double z1, double b );
double relbloch( double z12, double b1, double lambda, double theta0 );
double lindhard( double zz, double aa, double bb );
double Fbrems( double x );
double range( double e, double z1, double a1, int tno );
double qrange( double e, double z1, double a1, int tno );
double benton( double e, double z1, double a1, int tno );
double renergy( double e, double r0, double z1, double a1, int tno );
void run_range( FILE *finput, FILE *foutput );
void init_tables( int *initstat );
