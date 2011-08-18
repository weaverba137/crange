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
#include <math.h>

/*
 * Header file for complex arithmetic.
 * 
 */

typedef struct FCOMPLEX {double r, i;} fcomplex;

/* 
 * Declare functions from complex.c
 */
fcomplex Cadd( fcomplex a, fcomplex b );  /* Add two complex numbers */
fcomplex Csub( fcomplex a, fcomplex b );  /* Subtract two complex numbers */
fcomplex Cmul( fcomplex a, fcomplex b );  /* Multiply two complex numbers */
fcomplex Cdiv( fcomplex a, fcomplex b );  /* Divide two complex numbers */
fcomplex Complex( double re, double im ); /* Return complex given real 
                                             and imaginary parts */
fcomplex Conjg( fcomplex z );             /* Complex conjugation */
fcomplex Csqrt( fcomplex z );             /* Square root of a complex number */
fcomplex RCmul( double x, fcomplex a );   /* Multiply complex and 
                                             real numbers */
fcomplex Cexp( fcomplex z );              /* Complex exponential function */
fcomplex Cln( fcomplex z );               /* Complex logarithmic function */

fcomplex Cpow( fcomplex z, fcomplex c );  /* Complex number to complex power */
fcomplex Csin( fcomplex z );              /* Complex sine function */
fcomplex Ccos( fcomplex z );              /* Complex cosine function */
fcomplex Csinh( fcomplex z );             /* Complex sinh function */
fcomplex Ccosh( fcomplex z );             /* Complex cosh function */
double Cabs( fcomplex z );                /* Modulus of a complex number */
double Carg( fcomplex z );                /* Arguement of a complex number 
					     (in radians) */
double CRe( fcomplex z );                 /* Real part */
double CIm( fcomplex z );                 /* Imaginary part */
fcomplex Chyperg( fcomplex a, fcomplex b, fcomplex z );
                                          /* Complex confluent 
					     hypergeometric function */
fcomplex Clngamma( fcomplex z );          /* Complex logarithm of 
					     Complex Gamma function */
void multi_lngamma(int *NN, double *sN, double *eN, double *labsN, 
		   double *largN);        /* Utility function for 
					     multiple calls to Clngamma() */
void multi_hyperg( int *NN, double *arN, double *aiN, 
		   double *brN, double *biN, double *zrN, double *ziN, 
		   double *hyperrN, double *hyperiN );   
                                          /* Utility function for 
					     multiple calls to Chyperg() */

