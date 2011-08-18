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
 * Substantial portions of this file are based on Numerical Recipes in C,
 * 2nd ed., by W. H. Press, S. A. Teukolsky, W. T. Vetterling, and
 * B. P. Flannery, Cambridge: Cambridge University Press, 1992.
 *
 */
#include <math.h>
#include "complex.h"

/*
 * Miscellaneous defines
 */
#ifdef PI
#undef PI
#endif
#ifndef PI
#define PI 3.14159265359
#endif

/*
 * This package performs complex arithmetic, and has functions for the
 * fully complex Gamma function and confluent hypergeometric function. 
 * Many of these functions are found in W. H. Press, S. A. Teukolsky,
 * W. T. Vetterling, B. P. Flannery, Numerical Recipes in C, 2nd ed.,
 * (Cambridge University  Press, Cambridge, 1992).  In addition a number 
 * of useful analytic complex functions such as the exponential, logarithm,
 * trigonometric, and hyperbolic functions, are defined here.  Some further
 * description can be found in the header file "complex.h".  The branch
 * cut for the logarithm function is on the negative real axis. 
 */

fcomplex Cadd( fcomplex a, fcomplex b )
{
  fcomplex c;
  c.r=a.r+b.r;
  c.i=a.i+b.i;
  return c;
}

fcomplex Csub( fcomplex a, fcomplex b )
{
  fcomplex c;
  c.r=a.r-b.r;
  c.i=a.i-b.i;
  return c;
}

fcomplex Cmul( fcomplex a, fcomplex b )
{
  fcomplex c;
  c.r=a.r*b.r-a.i*b.i;
  c.i=a.i*b.r+a.r*b.i;
  return c;
}

fcomplex Cdiv( fcomplex a, fcomplex b )
{
  fcomplex c;
  double r, den;
  if (fabs(b.r) >= fabs(b.i)) {
    r=b.i/b.r;
    den=b.r+r*b.i;
    c.r=(a.r+r*a.i)/den;
    c.i=(a.i-r*a.r)/den;
  } else {
    r=b.r/b.i;
    den=b.i+r*b.r;
    c.r=(a.r*r+a.i)/den;
    c.i=(a.i*r-a.r)/den;
  }
  return c;
}

fcomplex Complex( double re, double im )
{
  fcomplex c;
  c.r=re;
  c.i=im;
  return c;
}

fcomplex Conjg( fcomplex z )
{
  fcomplex c;
  c.r=z.r;
  c.i= -z.i;
  return c;
}

fcomplex Csqrt( fcomplex z )
{
  fcomplex c;
  double x,y,w,r;
  if ((z.r == 0.0) && (z.i == 0.0)) {
    c.r=c.i=0.0;
    return c;
  } else {
    x=fabs(z.r);
    y=fabs(z.i);
    if (x >= y) {
      r=y/x;
      w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
    } else {
      r=x/y;
      w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
    }
    if (z.r >= 0.0) {
      c.r=w;
      c.i=z.i/(2.0*w);
    } else {
      c.i=(z.i >= 0.0) ? w : -w;
      c.r=z.i/(2.0*c.i);
    }
    return c;
  }
}

fcomplex RCmul( double x, fcomplex a )
{
  fcomplex c;
  c.r=x*a.r;
  c.i=x*a.i;
  return c;
}

fcomplex Cexp( fcomplex z )
{
  fcomplex ans;

  ans.r=exp(z.r)*cos(z.i);
  ans.i=exp(z.r)*sin(z.i);
  return ans;
}

fcomplex Cln( fcomplex z )
{
  fcomplex ans;

  ans.r=log(Cabs(z));
  ans.i=Carg(z);
  return ans;
}

fcomplex Cpow( fcomplex z, fcomplex c )
{
  fcomplex ans;

  ans = Cexp(Cmul(c,Cln(z)));
  return ans;
}

fcomplex Csin( fcomplex z )
{
  fcomplex iz, miz;

  iz = Cmul(Complex(0.0,1.0),z);
  miz = Cmul(Complex(0.0,-1.0),z);
  return(Cdiv(Csub(Cexp(iz),Cexp(miz)),Complex(0.0,2.0)));
}

fcomplex Ccos( fcomplex z )
{
  fcomplex iz, miz;

  iz = Cmul(Complex(0.0,1.0),z);
  miz = Cmul(Complex(0.0,-1.0),z);
  return(Cdiv(Cadd(Cexp(iz),Cexp(miz)),Complex(2.0,0.0)));
}

fcomplex Csinh( fcomplex z )
{
  fcomplex mz;

  mz=Cmul(Complex(-1.0,0.0),z);
  return(Cdiv(Csub(Cexp(z),Cexp(mz)),Complex(2.0,0.0)));
}

fcomplex Ccosh( fcomplex z )
{
  fcomplex mz;

  mz=Cmul(Complex(-1.0,0.0),z);
  return(Cdiv(Cadd(Cexp(z),Cexp(mz)),Complex(2.0,0.0)));
}

double Cabs( fcomplex z )
{
  double x,y,ans,temp;
  
  x=fabs(z.r);
  y=fabs(z.i);
  if (x == 0.0)
    ans=y;
  else if (y == 0.0)
    ans=x;
  else if (x > y) {
    temp=y/x;
    ans=x*sqrt(1.0+temp*temp);
  } else {
    temp=x/y;
    ans=y*sqrt(1.0+temp*temp);
  }
  return ans;
}

double Carg( fcomplex z )
{
  if(z.i==0.0){
    if(z.r >=0){
      return(0.0);
    } else {
      return(PI);
    }
  } else if(z.r==0.0) {
    if(z.i==0){
      return(0.0);
    } else if (z.i > 0.0) {
      return(PI/2.0);
    } else {
      return(-PI/2.0);
    }
  } else {
    return(atan2(z.i,z.r));
  }
}

double CRe( fcomplex z )
{
  return z.r;
}

double CIm( fcomplex z )
{
  return z.i;
}

fcomplex Chyperg( fcomplex a, fcomplex b, fcomplex z )
/* 
 * Confluent hypergeometric function.  WARNING, may not be stable for 
 * large values of |z|.
 */
{
  fcomplex Cone,Cm;
  fcomplex previousterm, term,sumterm;
  double dm;
  int m;

  Cone=Complex(1.0,0.0);
  term=Cone;
  sumterm=Cone;
  m=0;
  do {
    previousterm=term; 
    m+=1;
    dm=(double)m;
    Cm=Complex(dm-1.0,0.0);
    term=Cmul(previousterm,Cmul(Cdiv(Cadd(a,Cm),Cadd(b,Cm)),RCmul(1.0/dm,z)));
    sumterm=Cadd(sumterm,term);
  } while( Cabs(term) > 1.0e-6 && Cabs(previousterm) > 1.0e-6 );
  return(sumterm);
}

fcomplex Clngamma( fcomplex z )
/*
 * Fully complex logarithm of fully complex Gamma function.  Functional 
 * in all portions of the complex plane, including the negative real axis.
 * Note however, that the Gamma function has poles at all integers <= 0.
 */
{
  fcomplex result;
  double x, y, r, fj, cterm;
  double aterm1,aterm2,aterm3;
  double lterm1,lterm2,lterm3;
  int j;
  double num,denom;
  static double coeff[6]={76.18009172947146,
			  -86.50532032941677, 
			  24.01409824083091, 
			  -1.231739572450155, 
			  0.1208650973866179e-2, 
			  -0.5395239384953e-5};
  
  if(z.r>0) {
    x=z.r-1.0;
    y=z.i;
  } else {
    x=-z.r;
    y=-z.i;
  }
  r=sqrt((x+5.5)*(x+5.5)+y*y);
  aterm1=y*log(r);
  aterm2=(x+0.5)*atan2(y,(x+5.5))-y;
  lterm1=(x+0.5)*log(r);
  lterm2=-y*atan2(y,(x+5.5)) - (x+5.5) + 0.5*log(2.0*PI);
  num=0.0;
  denom=1.000000000190015;
  for(j=1;j<7;j++){
    fj=(double)j;
    cterm=coeff[j-1]/((x+fj)*(x+fj)+y*y);
    num+=cterm;
    denom+=(x+fj)*cterm;
  }
  num*=-y;
  aterm3=atan2(num,denom);
  lterm3 = 0.5*log(num*num + denom*denom);
  result.r=lterm1+lterm2+lterm3;
  result.i=aterm1+aterm2+aterm3;
  if(z.r<0){
    result=Csub(Complex(log(PI),0.0),Cadd(result,Cln(Csin(RCmul(PI,z)))));
  }
  return(result);
}

void multi_lngamma(int *NN, double *sN, double *eN, double *labsN, 
		   double *largN)
{
  fcomplex z,Clg;
  int j;
  double *labsout,*argout;

  for( j=0; j<*NN; j++){
    z=Complex( *(sN+j), *(eN+j) );
    labsout=labsN+j;
    argout=largN+j;
    Clg=Clngamma(z);
    *labsout=Clg.r;
    *argout=Clg.i;
  }
  return;
}

void multi_hyperg( int *NN, double *arN, double *aiN, 
		   double *brN, double *biN, double *zrN, double *ziN, 
		   double *hyperrN, double *hyperiN )
{
  fcomplex a,b,z,out;
  int j;
  double *hyperrout, *hyperiout;

  for( j=0; j<*NN; j++){
    a=Complex( *(arN+j), *(aiN+j) );
    b=Complex( *(brN+j), *(biN+j) );
    z=Complex( *(zrN+j), *(ziN+j) );
    hyperrout=hyperrN+j;
    hyperiout=hyperiN+j;
    out=Chyperg(a,b,z);
    *hyperrout=out.r;
    *hyperiout=out.i;
  }
  return;
}
