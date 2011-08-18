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
 * Include headers for variable and function declarations, &c..
 */
#include <crange.h>

int main( int argc, char **argv )
{
    FILE *finput,*foutput;
    short sswitch;
    int ini=0, have_switch=0, have_target=0;
    char inputname[50],outputname[50];
    char *switchfile, *targetfile;
    int c, errflag=0;
    /*
     * External variables used by getopt()
     */
    extern char *optarg;
    extern int optind, optopt;

    while((c=getopt(argc,argv,":hs:t:")) != -1) {
        switch (c) {
        case 'h':
            /*
             * Print help and exit.
             */
            errflag++;
            break;
        case 's':
            /*
             * Use filename as a switch file.
             */
            switchfile = optarg;
            have_switch++;
            break;
        case 't':
            /*
             * Use filename as a target file.
             */
            targetfile = optarg;
            have_target++;
            break;
        case ':':
            fprintf(stderr,"Option -%c requires an operand\n", optopt);
            errflag++;
            break;
        case '?':
            fprintf(stderr,"Unrecognised option: -%c\n", optopt);
            errflag++;
        }
    }
    if (errflag) {
        fprintf(stderr,"usage: crange [-h] [-s switch.ini] [-t target.dat] input [output]\n");
        return(1);
    }
    sswitch = (have_switch) ? init_switch(switchfile) : FD | FNS;
    ini = init_tables(targetfile);
    if(ini){
        fprintf(stderr,"Error initializing crange!\n");
        return(1);
    }
    if(argc-optind >= 1) {
        sscanf(argv[optind],"%s",inputname);
        finput=fopen(inputname, "r");
        if (finput==NULL) {
            fprintf(stderr,"Error opening task file!\n");
            return(2);
        }
    } else {
        fprintf(stderr,"No task file specified!\n");
        return(2);
    }
    if(argc-optind == 2) {
        sscanf(argv[optind+1],"%s",outputname);
        foutput=fopen(outputname, "w");
        if (foutput==NULL) {
            fprintf(stderr,"Error opening output file!\n");
            return(4);
        }
    } else {
        foutput=stdout;
    }
    run_range( finput, foutput, sswitch );
    fclose(finput);
    fclose(foutput);
    return(0);
}

double dedx( double e1, double rel0, double z0, double a1, short sswitch, int tno )
/*
 * This is the core of the whole package, the dE/dx calculator.  I have
 * based this largely on the work of M. H. Salamon, L.B.L. Report #10446
 * (1980).  Values of certain physical constants have been updated,
 * as well as some of the corrections to the basic stopping power formula.
 * Output is dE/dx in units of A MeV g^-1 cm^2.
 *
 * e1:   kinetic energy in A MeV
 * rel0: restricted energy loss parameter in eV
 * z0:   bare projectile charge
 * a1:   projectile atomic number
 * tno:  number of the current target (for looking up in the tdata struct)
 */
{
    extern tdata t[MAXAB];
    static double fva[10]={0.33,0.078,0.03,0.014,0.0084,
        0.0053,0.0035,0.0025,0.0019,0.0014};
    int i;
    double emass=0.511003e+6; /* eV/c^2 */
    double z2,a2;
    double g,b2,b,z1,z23;
    double delt;
    double etam2,cadj;
    double v,fv;
    double S,REL,f1,f2,f3,f4,f6,f7,f8,f9;
    double Sbr=0.0,Bbr;
    double Spa=0.0,dpa,ldpa,l0,Lpa0,Lpa0s,Lpa1,Lpa;
    double capA, capB;

    g=1.0+e1/931.4943;
    delt = ( sswitch & FD ) ? delta(g,tno) : olddelta(g,tno);
    b2=1.0-1.0/(g*g);
    b=sqrt(b2);
    z2=t[tno].z2;
    a2=t[tno].a2;
    /*
     * This is the modification of projectile charge due to electron
     * capture.  F. Hubert, R. Bimbot, and H. Gauvin, At. Data Nuc. Data
     * Tab. 46 (1990) 1, give an empirically determined function which
     * depends on the target material.  Two older versions, from
     * J. M. Anthony and W. A. Landford, Phys. Rev. A 25 (1982) 1868,
     * and T. E. Pierce and M. Blann, Phys. Rev. 173 (1968) 390 are also
     * available.
     */
    z23 = exp((2.0/3.0)*log(z0));
    if( sswitch & FE ) {
        if (z2 == 4.0) {
            z1=z0*(1.0 -
                ((2.045)+ 2.000*exp(-0.04369*z0))*
                exp(
                    -(7.000)*
                    exp((0.2643)*log(e1))*
                    exp(-(0.4171)*log(z0))
                    )
                );
        } else if (z2 == 6.0) {
            z1=z0*(1.0 -
                ((2.584)+ 1.910*exp(-0.03958*z0))*
                exp(
                    -(6.933)*
                    exp((0.2433)*log(e1))*
                    exp(-(0.3969)*log(z0))
                    )
                );
        } else {
            z1=z0*(1.0 -
                ((1.164 + 0.2319*exp(-0.004302*z2))+ 1.658*exp(-0.05170*z0))*
                exp(
                    -(8.144+0.9876*log(z2))*
                    exp((0.3140+0.01072*log(z2))*log(e1))*
                    exp(-(0.5218+0.02521*log(z2))*log(z0))
                    )
                );
        }
    } else {
        capA=1.16-z2*(1.91e-03 - 1.26e-05*z2);
        capB=(1.18-z2*(7.5e-03 - 4.53e-05*z2))/ALPHA;
        /*
         * capA=1.0;
         * capB=130.0;
         *
         * The Pierce and Blann formula can be activated by replacing the
         * variables capA and capB with the commented out values above.
         */
        z1=z0*(1.0-capA*exp(-capB*b/z23));
    }
    f1=0.3070722*z1*z1*z2/(b2*a1*a2);
    f2=log(2.0*emass*b2/t[tno].iadj);
    if( sswitch & FSH ){
        /*
         * The inner shell correction is somewhat problematic.  It arises when
         * the projectile velocity is comparable to the velocity of inner shell
         * electrons in the target medium.  This is discussed by U. Fano, Ann. Rev
         * Nucl. Sci. 13 (1963) 1. The shell correction can be included explicitly
         * using this formula from W. H. Barkas and M. J. Berger, Studies of
         * Penetration of Charged Particles in Matter, Natl. Acad. Sci. Pub. 1133
         * (1964).  Alternatively, the shell correction can be "hidden" in
         * the logarithmic mean ionization potential.  Much more work is required
         * before this topic can be fully understood.
         */
        etam2=1.0/(b*b*g*g);
        cadj=1.0e-6*(t[tno].iadj)*(t[tno].iadj)*etam2*(0.422377
            +etam2*(0.0304043-etam2*0.00038106))
            +1.0e-9*(t[tno].iadj)*(t[tno].iadj)*(t[tno].iadj)*etam2*(3.858019
            +etam2*(-0.1667989 + etam2*0.00157955));
        f2-=cadj/(z2);
        if( sswitch & FLE ){
            /*
             * The Leung, or relativistic shell correction is a small effect which
             * is due to relativistic inner shell electrons in very heavy targets.
             * See P. T. Leung, Phys. Rev. A 40 (1989) 5417, and P. T. Leung,
             * Phys. Rev. A 60 (1999) 2562.
             */
            f2-=(5.0/3.0)*log(2.0*emass*b2/t[tno].iadj)*(1.0e+03*t[tno].bind/(z2*emass))-(t[tno].iadj*t[tno].iadj/(4.0*emass*emass*b2));
        }
    }
    f6=2.0*log(g)-b2;
    /*
     * The Lindhard-Sorensen effect is now on by default.  The
     * Bloch-Mott-Ahlen effects are included for historical interest and
     * can be turned on by uncommenting the line after the next.
     */
    f3=lindhard(z1,a1,b,sswitch); /* comment out this line if uncommenting the next */
    /* f3=bma(z1,b); */
    f4=1.0;
    f8=0.0;
    f9=0.0;
    if( sswitch & FKIN ){
        /*
         * This an estimate of the ultrarelativistic kinematic correction from
         * S. P. Ahlen, Rev. Mod. Phys. 52 (1980) 121. It corrects to the
         * finite mass (as opposed to size) of the nucleus in relativistic
         * electron-nucleus collisions.
         */
        f8=0.5*(-log(1.0+2.0*((5.4858e-04)*g/a1)) - ((5.4858e-04)*g/a1)*b2/(g*g));
    }
    if( sswitch & FRAD ){
        /*
         * This is the radiative correction discussed in S. P. Ahlen,
         * Rev. Mod. Phys. 52 (1980) 121.  It arises from bremsstrahlung
         * of scattered electrons in ultrarelativistic collisions.  The
         * form here is that of V. Z. Jankus, Phys. Rev. 90 (1953) 4.
         * The parameter Q from that paper is here set equal to the geometric
         * mean between the the electron rest energy and 2 m_e c^2 gamma.
         */
        f9=(ALPHA/M_PI)*b2*(6.0822
            + log(2.0*g)*(
                log(2.0*g)*(
                    2.4167
                    + 0.3333*log(2.0*g))-8.0314));
    }
    if( sswitch & FPA ){
        /*
         * Slowing due to pair production.  This value and the value for
         * the bremsstrahlung correction below are based on the work of
         * A. H. Sorensen, "Stopping of relativistic heavy ions; the pair
         * production and bremsstrahlung channels", J. L. Duggan & I. L. Morgan
         * (Eds.), 17th Intl. Conf. on Application of Accelerators in Research
         * and Industry, Denton TX, Nov. 2002, AIP Conf. Proceedings No. 680,
         * pp 102-105.
         */
        dpa=1.0/sqrt(g);
        ldpa=log(dpa);
        l0=log(2.0*g);
        Lpa0=(19.0/9.0)*(log(g/4.0) - 11.0/6.0);
        Lpa0s=(19.0/9.0)*log(183.0*exp(-1.0/3.0*log(z2))/(1.0 + 4.0*6.25470095193633*183.0*exp(-1.0/3.0*log(z2))/g));
        Lpa1=dpa*(4178.0/(81*M_PI*M_PI) - 21.0/27.0 - 248.0*l0/(27.0*M_PI*M_PI)
            +(28.0*l0/9.0 - 446.0/27.0)*2.0*ldpa/(M_PI*M_PI) + 14.0*4.0*ldpa*ldpa/(9.0*M_PI*M_PI));
        Lpa=Lpa0s+Lpa1;
        Spa=4.08803936906434e-06*(z1*z1/a1)*(z2*z2/a2)*(1.0 + 1.0/z2)*g*Lpa;
    }
    if( sswitch & FBR ){
        /*
         * Slowing due to projectile bremsstrahlung.
         */
        Bbr=log(1.0 + 2.0*g*0.179524783764566/(exp((1.0/3.0)*log(a1)) + exp((1.0/3.0)*log(a2)))/a1);
        Sbr=5.21721169334564e-07*(z1*z1/a1)*(z1*z1/a1)*(z2*z2/a2)*g*Bbr;
    }
    if( sswitch & FBA ){
        /*
         * This is the Barkas correction as calculated in J. D. Jackson and
         * R. L. McCarthy, Phys. Rev. B 6 (1972) 4131.  It is multiplied
         * by a factor of two to bring it into agreement with J. Lindhard,
         * Nucl. Inst. Meth. 132 (1976) 1.  It is not, however, equal to the
         * results of Lindhard, and more work is needed to decide which, if any,
         * form is correct.  The recommended value seems to be the Jackson
         * & McCarthy result multiplied by two.  Jackson and McCarthy do not
         * have reliable values of F(V) for V < 0.8.  For the purposes of the
         * computation, the cut-off is placed at v=1.0.  I have followed the
         * convention of Salamon in having the Barkas correction multiply just
         * the "Bethe" portion of the stopping logarithm rather than the whole
         * stopping logarithm.  As there is considerable disagreement in the
         * literature about the application of correction, and as changing
         * the convention makes makes a difference of less than 1 A MeV even
         * in calculating the energy of stopping uranium, I have chosen to
         * leave it where it is.  Furthermore, I have found that a simple
         * power law v^-2 is adequate to model Jackson and McCarthy's function
         * for v > 1.0, so I have used this instead of the numbers found by
         * reading off one of Jackson and McCarthy's figures.
         */
        v=b*g/(ALPHA*sqrt(z2));
        if(v>1.0){
            /*
            if(v<9.0){
                i=1;
                while(v >= (double)i) i++;
                fv=fva[i-1]+(v-(double)(i-1))*(fva[i]-fva[i-1]);
            } else {
                fv=fva[9]*exp(-2.5*log(v/10.0));
            }
            */
            i=9;
            fv=fva[i]*exp(-2.0*log(v/10.0));
            f4=1.0 + 2.0*z1*fv/(sqrt(z2));
        } else {
            f4=1.0;
        }
        /*
        printf("%f %f\n",v,f4);
        */
    }
    S=f1*(f2*f4+f3+f6-(delt/2.0)+f8+f9) + Sbr + Spa;
    /*
     * Compute restricted energy loss.  REL is activated by setting rel0 >0.
     */
    if(rel0 > 0.0) {
        f7=log(2.0*emass*b2*g*g/rel0)+b2*(rel0/(2.0*emass*b2*g*g)-1.0);
        REL=f1*(f2*f4+f3+f6-delt/2.0 - 0.5*f7 +f8);
        return(REL);
    } else {
        return(S);
    }
}

double delta( double g, int tno )
/*
 * This function implements the density effect correction as formulated
 * in R. M. Sternheimer and R. F. Peierls, Phys. Rev. B 3 (1971) 3681 and as
 * extended in R. M. Sternheimer, M. J. Berger, S. M. Seltzer, Atom. Data
 * and Nucl. Data Tables 30 (1984) 261.  This version can distinguish
 * between solids and gasses, and between metals and insulators.  For
 * conducting materials, there is a low-energy density effect.
 */
{
    extern tdata t[MAXAB];
    double b,cbar,X,X0,X1;

    X0=t[tno].X0;
    X1=t[tno].X1;
    cbar=2.0*log(t[tno].iadj/t[tno].pla)+1.0;
    b=sqrt(1.0 - 1.0/(g*g));
    X=log10(b*g);
    if(t[tno].etad>0) {
        cbar-=2.303*log10(t[tno].etad);
        X1-=0.5*log10(t[tno].etad);
        X0-=0.5*log10(t[tno].etad);
    }
    if(X < X0) {
        return(t[tno].d0*exp(4.6052*(X-X0)));
    } else if (X >= X0 && X < X1) {
        return(4.6052*X + exp(log(t[tno].a) + t[tno].m*log(X1-X)) - cbar);
    } else {
        return( 4.6052*X - cbar );
    }
}

double olddelta( double g, int tno )
/*
 * This function implements the density effect correction as originally
 * formulated in R. M. Sternheimer and R. F. Peierls, Phys. Rev. B 3 (1971)
 * 3681.  Although it is now obsolete, I have included it here for
 * compatibility with earlier codes.
 */
{
    extern tdata t[MAXAB];
    double b,cbar,y;
    double y0,y1,dy3,a;

    if( g < 1.8 ) return(0.0);
    cbar=2.0*log(t[tno].iadj/t[tno].pla)+1.0;
    b=sqrt(1.0 - 1.0/(g*g));
    y=2.0*log(b*g);
    if( t[tno].etad > 0 ){
        y+=log(t[tno].etad);
        if(cbar>=12.25){
            y1=23.03;
            y0 = ( cbar >= 13.804 ) ? 1.502*cbar-11.52 : 9.212;
        } else {
            y1=18.42;
            if(cbar<12.25)y0=9.212;
            if(cbar<11.5)y0=8.751;
            if(cbar<11.0)y0=8.291;
            if(cbar<10.5)y0=7.830;
            if(cbar<10.0)y0=7.370;
        }
    } else {
        if(t[tno].iadj>=100.0){
            y1=13.82;
            y0 = ( cbar >= 5.215 ) ? 1.502*cbar-6.909 : 0.9212;
        } else {
            y1=9.212;
            y0 = ( cbar >= 3.681 ) ? 1.502*cbar-4.606 : 0.9212;
        }
    }
    if( y < y0 ) return( 0.0 );
    else {
        if( y > y1 ) return(y-cbar);
        else {
            dy3=(y1-y0)*(y1-y0)*(y1-y0);
            a=(cbar-y0)/dy3;
            return(y-cbar+a*(y1-y)*(y1-y)*(y1-y));
        }
    }
}

double bma( double z1, double b )
/*
 * This function is provided for historical reasons and implements the
 * Bloch, Mott & Ahlen corrections described below.
 */
{
    int n,msum;
    double b2,y,y2,sumr,fn,fn2,f3,f5;
    double lambda, theta0, cosx, st;
    fcomplex Cz1,Cz2;

    /*
     * This section includes the Mott correction of S. P. Ahlen, Phys. Rev.
     * A 17 (1978) 1236, the Bloch correction of F. Bloch, Ann. Phys.
     * (Leipzig) 16 (1933) 285, and the Ahlen correction of S. P. Ahlen,
     * Phys. Rev. A 25 (1982) 1856.  All three of these corrections are
     * rendered obsolete by the Lindhard-Sorensen correction, and are
     * included here for historical interest and comparison with older
     * calculations.  Note also that the Mott correction has a severely
     * limited range of validity, especially for high charges.  It's so
     * bad it can render the calculation not just inaccurate, but
     * unphysical (dE/dx < 0) below about 10 A MeV for uranium.  Ahlen
     * recommends turning the Mott correction off for Z/beta > 100.
     * Here for Z/beta > 100 the Mott correction is given the value at
     * Z/beta = 100. This prescription is given by C. J. Waddington,
     * P. S. Freier, and D. J. Fixsen, Phys. Rev. A 28 (1983) 464.
     */
    y=z1*ALPHA/b;
    y2=y*y;
    msum=(int)(10.0*y) + 1;
    sumr=0.0;
    for(n=1;n<msum;n++){
        fn=(double)n;
        fn2=fn*fn;
        sumr+=(1.0/(fn2+y2) - 1.0/fn2)/fn;
    }
    /*
     * Bloch and Ahlen corrections. The variables lambda and theta0 are
     * free parameters in the Ahlen correction.  Theta0 also appears in
     * the Mott correction.  Here I have used Ahlen's recommended values.
     * An alternative formula, theta0=sqrt(ALPHA/(b*g*lambda)), is
     * suggested by C. J. Waddington, P. S. Freier, and D. J. Fixsen,
     * Phys. Rev. A 28 (1983) 464.
     */
    lambda=1.0;
    theta0=0.1;
    f3=-y2*(1.202+sumr) + relbloch(z1,b,lambda,theta0);
    Cz1=Complex(0.5,-y);
    Cz2=Complex(1.0,y);
    cosx=cos(2.0*(CIm(Clngamma(Cz1))+CIm(Clngamma(Cz2))));
    st=sin(theta0/2.0);
    b2=b*b;
    /*
     * The Mott term
     */
    /*
    f5=0.5*z1a*(b*(1.725+(0.52-2.0*st*st)*M_PI*cosx)+
        z1a*(3.246-0.451*b2+z1a*(1.522*b+0.987/b+
        z1a*(4.569-0.494*b2-2.696/b2+
        z1a*(1.254*b+0.222/b
        -1.170/b2/b)))));
    */
    /*
    if( (exp(9.0*log(y))/6.0) > fabs(f5/(f2*f4+f3+f5+f6-delt/2.0)) )
    */
    if( y > 100.0*ALPHA ) y = 100.0*ALPHA;
    f5=0.5*b2*y*((1.725+(0.52-2.0*st)*M_PI*cosx)
        +y*((3.246-0.451*b2)
            +y*((0.987+1.522*b2)
                +y*((-2.696+b2*(4.569-0.494*b2))
                    +y*(-1.170+b2*(0.222+1.254*b2))
                    )
                )
            )
        );
    return( f3 + f5 );
}

double relbloch( double z12, double b1, double lambda, double theta0 )
/*
 * This is the relativistic Bloch (or Ahlen) correction of S. P. Ahlen,
 * Phys. Rev. A 25 (1982) 1856.  The evaluation of this correction has
 * been enormously simplified by the use of fully complex arithmetic.
 */
{
    fcomplex Cl1, Cl2;
    fcomplex Cf1,Cf2,Cf3,Cf4,Cf5,Cf6,Cf7;
    fcomplex Ci,Cnu,Cth,Cabgl;
    double ge = 0.5772157;
    double abgl;
    double g,nu;
    double sigma,es,bloch2;

    nu=z12*ALPHA/b1;
    g=1.0/sqrt(1.0-b1*b1);
    abgl=ALPHA/(b1*g*lambda);
    sigma=CIm(Clngamma(Complex(1.0,nu)));
    Ci=Complex(0.0,1.0);
    Cnu=Complex(1.0, 2.0*nu);
    Cth=Complex(theta0/2.0,0.0);
    Cabgl=Complex(abgl,0.0);
    Cf1=Cdiv(Complex(0.0,-1.0),Cnu);
    Cf2=Cpow(Cth,Cnu);
    Cf3=Cpow(Cabgl,Cnu);
    Cf4=Cexp(Complex(0.0,2.0*sigma));
    Cf5=Cdiv(Complex(1.0,0.0),Cnu);
    Cf6=Csub(Complex(log(theta0/2.0),0.0),Cf5);
    Cf7=Cadd(Complex(log(4.0/abgl)+ge-1.0,0.0),Cf5);
    Cl1=Cmul(Cf1,Csub(RCmul(2.0,Cf2),Cmul(Cf3,Cf4)));
    Cl2=Cmul(Cf1,Cadd(RCmul(2.0,Cmul(Cf2,Cf6)),Cmul(Cf3,Cmul(Cf4,Cf7))));
    es=(M_PI*nu)*exp(M_PI*nu)/sinh(M_PI*nu);
    bloch2=(M_PI/2.0)*b1*b1*nu*es*(4.0*nu*log(2.0)*CRe(Cl1)
        + (nu*M_PI-1.0)*CIm(Cl1)
        + 2.0*nu*CRe(Cl2));
    return(bloch2);
}

double lindhard( double zz, double aa, double bb, short sswitch )
/*
 * This is the Lindhard-Sorensen correction including finite nuclear
 * size effects as described in J. Lindhard and A. H. Sorensen,
 * Phys. Rev. A 53 (1996) 2443.  The defined variable FNS will turn off the
 * nuclear size effect if it is set to zero.  For values of the Lorentz
 * factor above 10/R, where R is the nuclear size divided by the electron
 * Compton wavelength, the correction is set to its asymptotic value
 * which is described by A. H. Sorensen,in Photonic, Electronic and Atomic
 * Collisions; Invited Papers XX Int. Conf. on the Physics of Electronic
 * and Atomic Collisions, eds. F. Aumayr and H. Winter (World Scientific,
 * Singapore, 1998) 475.  This also avoids some difficulties with the
 * evaluation of the confluent hypergeometric function (A. H. Sorensen,
 * private communication).
 */
{
    fcomplex Cexir, Cexis, Cedr, Ceds, Cske, Cmske, Cgrgs;
    fcomplex Clamr,Clams;
    fcomplex Caar, Cbbr, Caas, Cbbs, Czzr;
    fcomplex Cpi,Cone;
    int n,i,max;
    double compton=3.05573356675e-3; /* 1.18 fm / Compton wavelength */
    double a3,eta,gg,rho,prh;
    double a0,b0,a1,gz,z1,an,bn,anm1,bnm1,anp1,nn,asum,bsum,signk;
    double frgr, fsgs, grgs, figi, H;
    double k=0.0,k0,l,sk;
    double r,dk[3];
    double dkm1=0.0,dmk=0.0,sd2,sdm2,sdd;
    double term1,term2,term3,sumterm;
    double pct;
    double lls;

    Cpi=Complex(M_PI,0.0);
    Cone=Complex(1.0,0.0);
    a3=exp(log(aa)/3.0);
    eta=zz*ALPHA/bb;
    gg=1.0/sqrt(1.0 - bb*bb);
    rho=a3*compton;
    prh=bb*gg*rho;
    n = 1;
    sumterm = 0.0;
    term1 = 0.0;
    term3 = 1.0;
    term2 = 0.0;
    pct = 1.0;
    if ((gg < 10.0/rho) || !(sswitch & FNS)) {
/*      while(fabs(pct) > 0.01) {  */
        while ( n < 100 ) {
            k0 = (double)n;
            max = n == 1 ? 3 : 2;
            for(i=0;i<max;i++){
                if(i==0)k=k0;
                if(i==1)k= -k0-1.0;
                if(i==2)k= -k0;
                signk = k/fabs(k);
                sk = sqrt(k*k - ALPHA*ALPHA*zz*zz);
                l= (k>0) ? k : -k-1.0;
                Cske=Complex(sk+1.0,eta);
                Cexir=Csqrt(Cdiv(Complex(k,-eta/gg),Complex(sk,-eta)));
                Cedr=Cmul(Cexir,
                    Cexp(Complex(0.0,
                        -CIm(Clngamma(Cske))+(M_PI/2.0)*(l-sk))));
                if( sswitch & FNS ){
                    Cmske=Complex(-sk+1.0,eta);
                    Cexis=Csqrt(Cdiv(Complex(k,-eta/gg),Complex(-sk,-eta)));
                    Ceds=Cmul(Cexis,
                        Cexp(Complex(0.0,
                            -CIm(Clngamma(Cmske)) + (M_PI/2.0)*(l+sk))));
                    Caar=Cske;
                    Caas=Cmske;
                    Cbbr=Complex(2.0*sk + 1.0,0.0);
                    Cbbs=Complex(-2.0*sk + 1.0,0.0);
                    Czzr=Complex(0.0,2.0*prh);
                    Clamr=Cmul(Cmul(Cexir,Cexp(Complex(0.0,-prh))),
                        Chyperg(Caar,Cbbr,Czzr));
                    Clams=Cmul(Cmul(Cexis,Cexp(Complex(0.0,-prh))),
                        Chyperg(Caas,Cbbs,Czzr));
                    grgs=CIm(Clamr)/CIm(Clams);
                    Cgrgs=Clngamma(Cbbs);
                    grgs*=exp( CRe(Clngamma(Cske)) - CRe(Clngamma(Cmske))
                        -CRe(Clngamma(Cbbr)) + CRe(Cgrgs)
                        +2.0*sk*log(2.0*prh));
                    if(cos(CIm(Cgrgs))<1.0)grgs*= -1.0;
                    if(fabs(grgs) > 1.0e-9) {
                        frgr=sqrt((gg-1.0)/(gg+1.0))*CRe(Clamr)/CIm(Clamr);
                        fsgs=sqrt((gg-1.0)/(gg+1.0))*CRe(Clams)/CIm(Clams);
                        gz=-1.0*signk*(rho*gg+1.5*ALPHA*zz);
                        z1=-1.0*signk*zz;
                        b0=1.0;
                        a0=(1.0 + 2.0*fabs(k))*b0/(rho-gz);
                        a1=0.5*(gz+rho)*b0;
                        an=a1;
                        anm1=a0;
                        bnm1=b0;
                        asum=a0;
                        bsum=b0;
                        nn=1.0;
                        do{
                            bn=((rho-gz)*an + ALPHA*z1*anm1/2.0)/(2.0*nn+2.0*fabs(k)+1.0);
                            anp1=((gz+rho)*bn - ALPHA*z1*bnm1/2.0)/(2.0*nn + 2.0);
                            asum+=an;
                            bsum+=bn;
                            nn+=1.0;
                            anm1=an;
                            an=anp1;
                            bnm1=bn;
                        } while(fabs(anm1/asum) > 1e-6 && fabs(bnm1/bsum) > 1e-6 );
                        figi= (k>0) ? asum/bsum : bsum/asum;
                        H=((frgr-figi)/(figi-fsgs))*grgs;
                    } else {
                        H=0.0;
                    }
                } else {
                    H=0.0;
                    Ceds=Complex(0.0,0.0);
                }
                r=1.0 + H*H + 2.0*H*(Cedr.r*Ceds.r + Cedr.i*Cedr.i);
                dk[i]=Carg(Cadd(Cedr,RCmul(H,Ceds)));
            }
            if(n>1)dk[2]=dmk;
            sdm2=sin(dk[2]-dk[1]);
            term1 = k0*(k0+1.0)*sdm2*sdm2/(eta*eta*(2.0*k0 + 1.0));
            if(n>1) {
                sd2=sin(dk[0]-dkm1);
                term1+=k0*(k0-1.0)*sd2*sd2/(eta*eta*(2.0*k0 - 1.0));
            }
            sdd=sin(dk[0]-dk[2]);
            term2 = k0*sdd*sdd/(eta*eta*(4.0*k0*k0 - 1.0));
            term3 = term1 - 1.0/k0;
            sumterm += term2 + term3;
            n += 1;
            dkm1=dk[0];
            dmk=dk[1];
            pct = (term2 + term3)/sumterm;
        }
    } else {
        sumterm = -log(prh) - 0.2; /* Asymptotic value of the LS correction */
    }
    lls = sumterm + 0.5*bb*bb;
    return(lls);
}

double Fbrems( double x )
{
    int n=1;
    double t=1.0,s=0.0;

    if (x == 1.0) {
        return(M_PI*M_PI/12.0);
    } else if (x < 1.0) {
        while ((fabs(t) > 0.0001) || (n < 10)) {
            t*=-x;
            s+=t/((double)(n*n));
            n++;
        }
        return(-s);
    } else {
        while ((fabs(t) > 0.0001) || (n < 10)) {
            t*=-1.0/x;
            s+=t/((double)(n*n));
            n++;
        }
        return(M_PI*M_PI/12.0 + 0.5*log(x)*log(x) + s);
    }
}

double range( double e, double z1, double a1, short sswitch, int tno )
/*
 * This function computes total range given initial energy.  The technique
 * is quite clever, in that if from one call to the next, the projectile
 * and target material parameters do not change, the calculation of
 * range is performed by table interpolation rather than direct integration.
 * The savings in calculation time can be enormous.  However, the
 * range of valid energies is limited by the size of the table.  The
 * function dE/dx is evaluated at most of the energies contained in the
 * array tenerg[].  Results are stored in the two dimensional array
 * trange[][].  The second index of trange[][] can have an arbitrary range
 * and should be set to whatever is most useful.  Certainly it should
 * be no smaller than the number of target materials being used.
 * Output is range in units of g cm^-2.  The arguments are the same as
 * those described for the dedx() function.
 */
{
    extern double tenerg[MAXE],trange[MAXE][MAXAB];
    static double z1p = 0.0, a1p = 0.0;
    int table=1,i;
    double de2,dr,dedx1,dedx2,dedx3,dedx4,e1,e2,e3,e4;
    double rel = 0.0;

    if( (z1==z1p) && (a1==a1p) ) {
        if( (trange[ MAXE - 1 ][tno] > 0) ) table = 0;
    } else {
        z1p=z1;
        a1p=a1;
    }
    if(table){
        i=0;
        while(tenerg[i] < 8.0) {
            e1=tenerg[i];
            trange[i][tno]=benton(e1, z1, a1, tno);
            i++;
        }
        while(i<MAXE){
            de2=(tenerg[i]-tenerg[i-1])/2.0;
            e1=tenerg[i-1]+1.33998104*de2;
            dedx1=dedx(e1,rel,z1,a1,sswitch,tno);
            e2=tenerg[i-1]+1.86113631*de2;
            dedx2=dedx(e2,rel,z1,a1,sswitch,tno);
            e3=tenerg[i-1]+0.13886369*de2;
            dedx3=dedx(e3,rel,z1,a1,sswitch,tno);
            e4=tenerg[i-1]+0.66001869*de2;
            dedx4=dedx(e4,rel,z1,a1,sswitch,tno);
            dr=de2*(0.65214515/dedx1 + 0.34785485/dedx2 + 0.34785485/dedx3
                + 0.65214515/dedx4);
            trange[i][tno] = trange[i-1][tno] + dr;
            i++;
        }
    }
    if( e > tenerg[0]) {
        i=1;
        while( e > tenerg[i] ) i++;
        return( trange[i-1][tno] + ( e - tenerg[i-1] )
            *(trange[i][tno]-trange[i-1][tno])/(tenerg[i]-tenerg[i-1]) );
    } else {
        return( e*trange[0][tno]/tenerg[0] );
    }
}

double qrange( double e, double z1, double a1, short sswitch, int tno )
/*
 * This function computes total range given initial energy.  This
 * is a direct integration of dE/dx from 1 A MeV up to the initial energy.
 */
{
    int i;
    double de2,dr,dedx1,dedx2,dedx3,dedx4,e1,e2,e3,e4;
    double ei=8.0,lef,lei,en=1.0,pen,rel = 0.0;
    double ln10,entries;
    double ra[MAXE];

    if(e > ei){
        ra[0]=benton(ei,z1,a1,tno);
        i=1;
        ln10=log(10.0);
        lei=log10(ei);
        lef=log10(e);
        entries = (double) (MAXE - 1);
        do {
            pen=en;
            en=exp(ln10*(lei + ((double)i)*lef/entries));
            de2=(en-pen)/2.0;
            e1=pen+1.33998104*de2;
            dedx1=dedx(e1,rel,z1,a1,sswitch,tno);
            e2=pen+1.86113631*de2;
            dedx2=dedx(e2,rel,z1,a1,sswitch,tno);
            e3=pen+0.13886369*de2;
            dedx3=dedx(e3,rel,z1,a1,sswitch,tno);
            e4=pen+0.66001869*de2;
            dedx4=dedx(e4,rel,z1,a1,sswitch,tno);
            dr=de2*(0.65214515/dedx1 + 0.34785485/dedx2 + 0.34785485/dedx3
                + 0.65214515/dedx4);
            ra[i] = ra[i-1] + dr;
            i++;
        }while(en < e);
        return(ra[i-1]);
    } else if (e > 1.0 && e <= ei) {
        return( benton(e,z1,a1,tno) );
    } else {
        return( e*benton(en,z1,a1,tno)/en );
    }
}

double benton( double e, double z1, double a1, int tno )
/*
 * This function is the result of  empirical fits to very low energy
 * 1 A MeV < E < 8 A MeV ion ranges.  It follows the methods of
 * W. H. Barkas and M. J. Berger, Studies of Penetration of Charged
 * Particles in Matter, Natl. Acad. Sci. Pub. 1133 (1964).  A simplified
 * discussion, with a more complicated formula is given in
 * E. V. Benton and R. P. Henke, Nucl. Inst. Meth. 67 (1969) 87.
 * As yet I know of no nicer way to deal with these low energies.
 */
{
    extern tdata t[MAXAB];
    int l,m,n;
    double g,b,bzz,x;
    double term,logt,logi,loglambda;
    double cr = 938.2723/931.4943;
    double prnglo[3];
    double cz[4],join[4];
    static double amn[3][4][4] = {
        {
            {-8.72500,  1.88000,  0.741900,  0.752000},
            { 0.83090,  0.11140, -0.528800, -0.555890},
            {-0.13396, -0.06481,  0.126423,  0.128431},
            { 0.01262,  0.00540, -0.009341, -0.009306}
        },
        {
            {-7.6604e-01,  2.5398e+00, -2.4598e-01,  0.0000e+00},
            { 7.3736e-02, -3.1200e-01,  1.1548e-01,  0.0000e+00},
            { 4.0556e-02,  1.8664e-02, -9.9661e-03,  0.0000e+00},
            { 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00}
        },
        {
            {-8.0155e+00,  1.8371e+00,  4.5233e-02, -5.9898e-03},
            { 3.6916e-01, -1.4520e-02, -9.5873e-04, -5.2315e-04},
            {-1.4307e-02, -3.0142e-02,  7.1303e-03, -3.3802e-04},
            { 3.4718e-03,  2.3603e-03, -6.8538e-04,  3.9405e-05}
        }
    };
    static double czmn[4][4] = {
        {-6.000e-05,  5.252e-02,  1.285e-01,  0.000e+00},
        {-1.850e-03,  7.355e-02,  7.171e-02, -2.723e-02},
        {-7.930e-02,  3.323e-01, -1.234e-01,  1.530e-02},
        { 2.200e-01,  0.000e+00,  0.000e+00,  0.000e+00}
    };
    static double cjoin[2][7] = {
        { 0.94,  20.19, -84.08,  132.98, -30.77, -102.29, 64.03},
        {12.62, -51.96, 199.14, -367.09, 327.06, -108.57,  0.00}
    };

    logt=log( e * cr );
    logi=log( t[tno].iadj );
    for (l = 0; l < 2; l++) {
        join[l] = cjoin[l][m=6];
        while (m > 0) join[l] = 0.001*t[tno].iadj*join[l] + cjoin[l][--m];
    }
    /*
     * join[l] demarcates the three energy regions represented by the
     * three sets of coefficients in amn.  The demarcation is variable
     * in order to minimize discontinuities at the boundary.  The
     * coefficients in cjoin are inhereted from legacy code; I have not found
     * them in the non-obscure literature.  Approximately, the three
     * regions are E < 1 A MeV, 1 < E < 7 A MeV, and E > 7 A MeV.
     */
    for (l = 0; l < 3; l++) {
        loglambda=0.0;
        for (m = 3; m >= 0; m-- ) {
            term=amn[l][m][n=3];
            while (n > 0) term=term*logt + amn[l][m][--n];
            loglambda+=term;
            loglambda*=logi;
        }
        loglambda/=logi;
        loglambda+=log( (t[tno].a2)/(t[tno].z2) );
        prnglo[l]=exp(loglambda);
        if (l == 1) prnglo[l]*=1.0e-03;
    }
    g=1.0 + e/931.4943;
    b=sqrt(1.0 -1.0/(g*g));
    x=137.0*b/z1;
    for (m = 0; m < 4; m++) {
        cz[m]=0.0;
        for (n = 3; n >= 0; n--) {
            cz[m]+=czmn[m][n];
            cz[m]*=x;
        }
        cz[m]/=x;
    }
    l = 1;
    if ( e < join[0]) l = 0;
    if ( e > join[1]) l = 2;
    if (x <= 0.2) {
        n = 0;
    } else if ( x > 0.2 && x <= 2.0 ) {
        n = 1;
    } else if ( x > 2.0 && x <= 3.0 ) {
        n = 2;
    } else if ( x > 3.0 ) {
        n = 3;
    }
    bzz=(31.8+3.86*exp((5.0/8.0)*logi))
        *( (t[tno].a2)/(t[tno].z2) )*1.0e-06*exp((8.0/3.0)*log( z1 ));
    return(( (a1/cr)/(z1*z1) )*(prnglo[l] + bzz*cz[n]));
}

double renergy( double e, double r0, double z1, double a1, short sswitch, int tno )
/*
 * This function finds extracts energies from the range tables by
 * table interpolation.  Inputs are range measured in g cm^-2 , and other
 * inputs described above. The first call to range() is to initialize
 * the range-energy tables or find the correct table if it has
 * already been computed.
 */
{
    extern double trange[MAXE][MAXAB],tenerg[MAXE];
    double rr,r;
    int i;

    if( e > 0.0 ){
        rr = range(e,z1,a1,sswitch,tno);
        r = rr - r0;
        if( r < 0.0 ) return(0.0);
    } else {
        rr = range(tenerg[0],z1,a1,sswitch,tno);
        r = r0;
    }
    if(r > trange[0][tno]){
        i=1;
        while( r > trange[i][tno] ) i++;
        return(tenerg[i-1]+(r-trange[i-1][tno])*(tenerg[i]-tenerg[i-1])/
            (trange[i][tno]-trange[i-1][tno]));
    } else {
        return(tenerg[0]*r/trange[0][tno]);
    }
}

void run_range( FILE *finput, FILE *foutput, short sswitch )
/*
 * The following is a utility function which steps through the
 * range, energy and dE/dx tasks specified in the input data file.
 * The tasks are denoted by a single letter:  r to compute ranges,
 * e to compute energies, and d to compute dE/dx.  The task letter should
 * be followed by the energy (or range) at which to compute range (or
 * energy), the charge and mass of the particle, and the name of the
 * target material.  Names of target materials can be found in the
 * target.dat file.  Target material names may be up to eight (8)
 * characters in length and should contain no whitespace.
 */
{
    extern tdata t[MAXAB];
    char task[2];
    char abs[9],previous[9];
    double red1,red2,z1,a1;
    double out=0.0;
    int icols=6;
    int k=0;

    for(;;){
        icols=fscanf(finput,"%s %lf %lf %lf %lf %s\n",
            task,&red1,&red2,&z1,&a1,abs);
        if(icols < 6) break;
        if( strcmp( abs, previous ) != 0 ) {
            k=0;
            while( (strcmp( t[k].tname, abs ) != 0)
                && (strcmp( t[k].tname, "Unknown" ) != 0)) k++;
        }
        strcpy( previous, abs );
        out=0.0;
        if(icols==6 && strcmp( t[k].tname, "Unknown" ) != 0){
            if(strcmp( task, "r" )==0){
                out=range(red1,z1,a1,sswitch,k);
            } else if(strcmp( task, "e" )==0){
                out=renergy(red1,red2,z1,a1,sswitch,k);
            } else if(strcmp( task, "d" )==0){
                out=dedx(red1,red2,z1,a1,sswitch,k);
            }
        }
        fprintf(foutput,"%f\n",out);
    }
    return;
}

short init_switch( char *switchfile )
/*
 * Initializes the value of sswitch.
 */
{
    short sswitch;
#ifdef HAVE_INIPARSER_H
    dictionary *d;
#endif
    sswitch = FD | FNS;
#ifdef HAVE_INIPARSER_H
    if (access(switchfile,R_OK)) {
        sswitch = 0;
        d = iniparser_load( switchfile );
        if (iniparser_getboolean(d,"crange:barkas",0)) sswitch |= FBA;
        if (iniparser_getboolean(d,"crange:shell",0))  sswitch |= FSH;
        if (iniparser_getboolean(d,"crange:leung",0))  sswitch |= FLE;
        if (iniparser_getboolean(d,"crange:new delta",1))  sswitch |= FD; /* True by default! */
        if (iniparser_getboolean(d,"crange:new electron capture",0))  sswitch |= FE;
        if (iniparser_getboolean(d,"crange:finite nuclear size",1))  sswitch |= FNS; /* True by default! */
        if (iniparser_getboolean(d,"crange:kinematic",0))  sswitch |= FKIN;
        if (iniparser_getboolean(d,"crange:radiative",0))  sswitch |= FRAD;
        if (iniparser_getboolean(d,"crange:pair",0))  sswitch |= FPA;
        if (iniparser_getboolean(d,"crange:bremsstrahlung",0))  sswitch |= FBR;
        iniparser_freedict(d);
    } else {
        fprintf(stderr,"Could not read switch file: %s. Using default crange switches.\n",switchfile);
    }
#endif
    return sswitch;
}

int init_tables( char* foo )
/*
 * This utility creates and sets to zero all entries in the external
 * absorber and range tables.  It also sets up the energy table by
 * creating a logarithmically uniform distribution of energies between
 * some minimum energy and some maximum energy, with a number of entries
 * given by MAXE.  Finally, it opens the absorber data file for read-only
 * access. The variable pointed to by initstat will be set to zero on
 * successful completion of the initialization.
 */
{
    extern tdata t[MAXAB];
    extern double trange[MAXE][MAXAB], tenerg[MAXE];
    int i,j,k,tcol,initstat;
    double ln10,l10Emin,l10Emax,decades,entries;
    char rswitch[20];
    char root1[50] = DEDX ;
    char *targetfile;
    FILE *fabsorber;

    targetfile = strcat( root1, "/target.dat" );
    ln10=log(10.0);
    l10Emin=0.0; /* minimum energy 1 A MeV */
    l10Emax=6.0; /* maximum energy 1 A TeV */
    decades=l10Emax-l10Emin;
    entries=MAXE - 1.0;
    for(i=0;i < MAXAB;i++){
        t[i].z2=t[i].a2=t[i].iadj=t[i].rho=t[i].pla=0.0;
        t[i].X0=t[i].X1=t[i].a=t[i].m=t[i].d0=t[i].etad=t[i].bind=0.0;
        for(j=0;j < MAXE;j++){
            if(i==0) tenerg[j]=exp(ln10*(l10Emin + ((double)j)*decades/entries));
            trange[j][i]=0.0;
        }
    }
    k=0;
    fabsorber=fopen( targetfile , "r" );
    do{
        tcol=fscanf(fabsorber,
            "%s %lf %lf %lf %le %le %lf %lf %lf %lf %lf %lf %le\n",
            t[k].tname,&(t[k].z2),&(t[k].a2),&(t[k].iadj),&(t[k].rho),
            &(t[k].pla),&(t[k].X0),&(t[k].X1),&(t[k].a),&(t[k].m),
            &(t[k].d0),&(t[k].etad),&(t[k].bind));
        if(tcol < 13) printf("Error at target %s :  %d items converted.\n",t[k].tname,tcol);
        k++;
    }while( strcmp( t[k-1].tname, "Unknown" ) != 0 );
    initstat=fclose(fabsorber);
    return initstat;
}

