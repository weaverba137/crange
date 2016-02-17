/**
 * \mainpage crange - The Berkeley Range-Energy Calculator
 * \author Benjamin Weaver <benjamin.weaver@nyu.edu>
 * \version 1.6.1
 *
 * \copyright (C) 2001-2011 Benjamin Weaver, LGPL
 *
 * \section sec-intro Introduction
 *
 * Thank you for choosing the Berkeley Range-Energy Calculator.
 *
 * \section sec-preinstall Pre-Installation
 *
 * \subsection subsec-sysreq System Requirements
 *
 * crange should run on any POSIX-based system that has the appropriate
 * libraries installed (see \ref subsec-prereq "Prerequisites").
 *
 * \subsection subsec-prereq Library Prerequisites
 *
 * \subsubsection subsubsec-gsl GNU Scientific Library
 *
 * The GNU Scientific Library (http://www.gnu.org/s/gsl/) is required.
 * This code was most recently tested with version 1.15, though any version
 * with full complex variable support should work.
 *
 * \subsubsection subsubsec-iniparser iniParser
 *
 * The iniParser Library (http://ndevilla.free.fr/iniparser/), version 3.0 or
 * later is strongly recommended.  The code will still compile if it is
 * not installed, but the resulting binary will not be able to read optional
 * target or switch information, so only compiled-in defaults will be
 * available.
 *
 * \section sec-install Installation
 *
 * See the INSTALL file.
 *
 * \section sec-running Running
 *
 * \subsection subsec-command On the Command Line
 *
 * Type crange -h to get the list of command-line options.
 *
 * \subsection subsec-tasks Task List
 *
 * The crange program executes a list of tasks from a file.  Here's an example:
 *
 * \verbinclude tasks.txt
 *
 * The first line tells the calculator to compute the range (in g cm<sup>-2</sup>) of
 * uranium (Z=92,A=238) at a kinetic energy of 1200 MeV per nucleon, in an
 * aluminum target.  The second line tells the calculator to compute the
 * initial kinetic energy (in A MeV) of gold (Z=79,A=197) whose range in
 * the plastic track-etch detector CR-39 was 9.2 g cm<sup>-2</sup>.  The third line asks
 * for the final kinetic energy after passing through 9.2 g cm<sup>-2</sup> of CR-39, given
 * an initial energy of 10.6 A GeV. The fourth line tells the calculator to
 * compute dE/dx (in A MeV g<sup>-1</sup> cm<sup>2</sup>) for gold with kinetic energy 10.6 A GeV in
 * air.  The fifth line computes REL instead of dE/dx with the REL cutoff set to
 * 300 eV. The sixth line computes primary ionization.
 * Arguments which are zero (0) are dummies which are necessary for
 * place-holding.
 *
 * The list of tasks may be of any length in any combination of ranges,
 * energies or dE/dx.  The order of the output will be the same as the
 * order of the input.
 *
 * \subsection subsec-switch Switching Optional Effects
 *
 * The switch.ini file included with the source distribution shows the
 * effects that may be turned on or off by the user.  The values in the
 * switch.ini file are the default values that are compiled into the program.
 * The user may supply a modified switch.ini file on the command line.
 * There are additional details in the file itself.
 *
 * \subsection subsec-target Adding or Modifying Targets
 *
 * The target.ini file included with the source distribution lists the
 * targets that are compiled in to the program by default.  Most of the target
 * data is taken from Sternheimer, Berger \& Seltzer, \cite art_rms1.
 * The definitions of the material properties are in the target.ini file.
 * The user may add additional targets or orverride existing target values
 * by supplying a different target.ini file (with the same format!) on
 * the command line.
 *
 * \section sec-updates Updates
 *
 * Visit http://sdss.physics.nyu.edu/~bw55/dedx/ .
 *
 * \section sec-history History
 *
 *  - 1.6.1: Fix version strings.
 *  - 1.6.0: Make crange compatible with GNU autotools and GNU Scientific Library.
 *  - 1.5.3: Simplified calculation switching and added a crange.h file
 *  - 1.5.2: Increased value of MAXAB to acommodate a larger target data file
 *  - 1.5.1: Fixed typo in one of the parameters of the electron capture corection
 *
 */
/**
 * \file crange.c
 * \brief Source code for crange.
 *
 * This file contains all source code for the crange executable.
 */
/*
 * Include headers for variable and function declarations, &c..
 */
#include <crange.h>
/**
 * \brief Main crange program.
 *
 * Standard C main() program.
 *
 * \param argc Number of command line options.
 * \param argv The command line options.
 *
 * \return The exit status.
 */
int main( int argc, char **argv )
{
    FILE *finput,*foutput;
    tdata *extratargets, *listdummy;
    short sswitch;
    int have_switch=0, have_target=0, have_command=0, have_output=0;
    char inputname[50];
    char *switchfile, *targetfile, *command, *outputname;
    int c, errflag=0, listflag=0, fd=-1;
    char tempfilename[15] = "";
    static char list[5] = "List";
    extern int errno; /* From errno.h */
    extern char *optarg; /* External variable used by getopt(). */
    extern int optind, optopt; /* External variable used by getopt(). */
    /* End declarations */
    while((c=getopt(argc,argv,":c:hlo:s:t:")) != -1) {
        switch (c) {
        case 'c':
            /*
             * Interpret the value of optarg as a command to the crange engine.
             */
             command = optarg;
             have_command++;
             break;
        case 'h':
            /*
             * Print help and exit.
             */
            errflag++;
            break;
        case 'l':
            /*
             * Print the built-in target table and exit.
             */
            listflag++;
            break;
        case 'o':
            /*
             * Print output to a file
             */
             outputname = optarg;
             have_output++;
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
        fprintf(stderr,"usage: crange [-c COMMAND] [-h] [-l] [-o FILE] [-s switch.ini] [-t target.ini] <task file>\n");
        fprintf(stderr,"       -c COMMAND    = Execute this one-line command instead of reading it from a file.\n");
        fprintf(stderr,"       -h            = Print this help message and exit.\n");
        fprintf(stderr,"       -l            = Print the built-in target table and exit.\n");
        fprintf(stderr,"       -o FILE       = Write to this file instead of standard output.\n");
        fprintf(stderr,"       -s switch.ini = Override the default switch values by reading this file.\n");
        fprintf(stderr,"       -t target.ini = Override the default target values by reading this file.\n");
        fprintf(stderr,"       <task file>   = A file containing a list of tasks for crange.  Required unless a command is specified with -c.\n");
        return(1);
    }
    if (listflag) {
        listdummy = find_target(list,NULL);
        return(0);
    }
    sswitch = (have_switch) ? init_switch(switchfile) : SSWITCH_DEFAULT;
    extratargets = (have_target) ? init_target(targetfile) : NULL;
    if(argc-optind >= 1) {
        sscanf(argv[optind],"%s",inputname);
        finput=fopen(inputname, "r");
        if (finput==NULL) {
            fprintf(stderr,"Error opening task file!\n");
            return(2);
        }
    } else if (have_command) {
        /*
         * Create a temporary file to hold the command.
         */
        strcpy(tempfilename, "/tmp/cr.XXXXXX");
        if ((fd = mkstemp(tempfilename)) == -1 || (finput = fdopen(fd, "w+")) == NULL) {
            if (fd != -1) {
                close(fd);
                unlink(tempfilename);
            }
            fprintf(stderr, "%s: %s\n", tempfilename, strerror(errno));
            return(2);
        }
        /*
         * Write the command to the temporary file.
         */
        fprintf(finput,"%s\n",command);
        rewind(finput);
    } else {
        fprintf(stderr,"No task file specified!\n");
        return(2);
    }
    if(have_output) {
        foutput=fopen(outputname, "w");
        if (foutput==NULL) {
            fprintf(stderr,"Error opening output file!\n");
            return(4);
        }
    } else {
        foutput=stdout;
    }
    init_table();
    run_range( finput, foutput, sswitch, extratargets );
    fclose(finput);
    fclose(foutput);
    if (have_command) unlink(tempfilename);
    if (have_target) free(extratargets);
    return(0);
}
/**
 * \brief Confluent hypergeometric function.
 *
 * Computes the confluent hypergeometric function.  All input parameters
 * are complex numbers.  Uses the formula:
 * \f[ M(a,b,z) = 1 + \sum_{n=1} \frac{(a)_n}{(b)_n}\frac{z^n}{n!} , \f]
 * where
 * \f[ (x)_n \equiv \frac{\Gamma(x+n)}{\Gamma(x)} \f]
 * is the Pochhammer Symbol.
 *
 * \param a First parameter of the hypergeometric function.
 * \param b Second parameter of the hypergeometric function.
 * \param z A complex number.
 *
 * \return The value \f$ M(a,b,z) \f$ , a complex number.
 *
 * \warning May not be stable for large values of \f$|z|\f$.
 */
gsl_complex complex_hyperg( gsl_complex a, gsl_complex b, gsl_complex z )
{
    gsl_complex Cm, previousterm, term, sumterm;
    double dm = 0.0;
    /* End declarations */
    term=GSL_COMPLEX_ONE;
    sumterm=GSL_COMPLEX_ONE;
    do {
        previousterm=term;
        dm+=1.0;
        Cm=gsl_complex_rect(dm-1.0,0.0);
        term=gsl_complex_mul(previousterm,
            gsl_complex_mul(
                gsl_complex_div(gsl_complex_add(a,Cm),
                    gsl_complex_add(b,Cm)),gsl_complex_div_real(z,dm)));
        sumterm=gsl_complex_add(sumterm,term);
    } while( gsl_complex_abs(term) > 1.0e-6 && gsl_complex_abs(previousterm) > 1.0e-6 );
    return(sumterm);
}
/**
 * \brief Complex logarithm of the Gamma function.
 *
 * Computes the fully complex logarithm of the fully complex Gamma function.
 * Works in all portions of the complex plane, including the negative real
 * axis.
 *
 * \param z A complex number.
 *
 * \return \f$ \ln \Gamma(z) \f$ , a complex number.
 *
 * \warning The Gamma function has poles at all integers \<= 0.
 */
gsl_complex complex_lngamma( gsl_complex z )
{
    gsl_complex result;
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
    /* End declarations */
    if(GSL_REAL(z)>0) {
        x=GSL_REAL(z)-1.0;
        y=GSL_IMAG(z);
    } else {
        x=-GSL_REAL(z);
        y=-GSL_IMAG(z);
    }
    r=sqrt((x+5.5)*(x+5.5)+y*y);
    aterm1=y*log(r);
    aterm2=(x+0.5)*atan2(y,(x+5.5))-y;
    lterm1=(x+0.5)*log(r);
    lterm2=-y*atan2(y,(x+5.5)) - (x+5.5) + 0.5*log(2.0*M_PI);
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
    GSL_SET_COMPLEX(&result,lterm1+lterm2+lterm3,aterm1+aterm2+aterm3);
    if(GSL_REAL(z)<0){
        result=gsl_complex_sub(gsl_complex_rect(log(M_PI),0.0),
            gsl_complex_add(result,
                gsl_complex_log(gsl_complex_sin(
                    gsl_complex_mul_real(z,M_PI)))));
    }
    return(result);
}
/**
 * \brief Computes effective projectile charge.
 *
 * This is the modification of projectile charge due to electron
 * capture.  Hubert, Bimbot \& Gauvin, \cite art_fh2, give an
 * empirically determined function which depends on the target material.
 * This version is used if #SSWITCH_EC is set.
 * Two older versions, from Anthony \& Landford, \cite art_jma,
 * and Pierce \& Blann, \cite art_tep are also available.
 *
 * \param z0 The bare projectile charge.
 * \param e1 The projectile kinetic energy in A MeV.
 * \param z2 The target mean nuclear charge.
 * \param sswitch The switch bit field.
 *
 * \return The effective projectile charge.
 *
 * \bug The Pierce \& Blann formula is not actually available; it is simply
 * commented out.
 */
double effective_charge( double z0, double e1, double z2, short sswitch )
{
    double z23, z1, g, b2, b;
    double capA, capB;
    /* End declarations */
    g=1.0+e1/ATOMICMASSUNIT;
    b2=1.0-1.0/(g*g);
    b=sqrt(b2);
    z23 = exp((2.0/3.0)*log(z0));
    if( sswitch & SSWITCH_EC ) {
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
    return(z1);
}
/**
 * \brief Computes primary ionization.
 *
 * This computes the primary ionization, the number of delta-rays produced
 * per unit length.  The formula is based on Bethe \cite art_hb,
 * as well as Fleischer <em>et al.</em>, \cite art_rlf3.
 *
 * \param e1 The projectile kinetic energy in A MeV.
 * \param z0 The projectile charge.
 * \param I0 The binding energy of outermost electron in eV.
 * \param f0 The fraction of electrons in the outermost state.
 * \param K  A constant that depends on the target.
 * \param sswitch The switch bit field.
 * \param target A pointer to a ::TDATA structure.
 *
 * \return Number of delta-rays per unit length in units of g<sup>-1</sup> cm<sup>2</sup>.
 *
 * \bug The parameters needed are not contained in the target table.
 */
double djdx( double e1, double z0, double I0, double f0, double K, short sswitch, tdata *target)
{
    double z2,a2;
    double g,b2,b,z1;
    double delt;
    double J,f1,f2;
    /* End declarations */
    g=1.0+e1/ATOMICMASSUNIT;
    delt = ( sswitch & SSWITCH_ND ) ? delta(g,target) : olddelta(g,target);
    b2=1.0-1.0/(g*g);
    b=sqrt(b2);
    z2=target->z2;
    a2=target->a2;
    z1 = effective_charge(z0, e1, z2, sswitch);
    f1=0.3070722*z1*z1*z2/(b2*a2);
    f2=log(2.0*ELECTRONMASS*b2*g*g/I0);
    J=0.5*f1*(f2-b2-delt+K)*(f0/I0);
    return(J);
}
/**
 * \brief Computes dE/dx.
 *
 * This is the core of the whole package, the dE/dx calculator.  I have
 * based this largely on the work of Salamon, \cite tech_mhs.
 * Values of certain physical constants have been updated,
 * as well as some of the corrections to the basic stopping power formula.
 *
 * If the restricted energy loss parameter \em rel0 is non-zero, dedx() computes
 * restricted energy loss instead.
 *
 * The dE/dx calculator includes a number of effects that are controlled by
 * switches encoded in a bit field.  Below we describe each bit field and the
 * effect it controls.
 *
 *  - #SSWITCH_ND : Density effect version. If this bit is set (which it is
 *    by default), a newer version of the density effect is used.  See
 *    delta() and olddelta() for details.
 *  - #SSWITCH_SH : Inner shell correction. The inner shell correction is
 *    somewhat problematic.  It arises when the projectile velocity is
 *    comparable to the velocity of inner shell electrons in the target medium.
 *    This is discussed by Fano, \cite art_uf. The shell correction can be
 *    included explicitly using this formula from Barkas \& Berger, \cite coll_whb.
 *    Alternatively, the shell correction can be "hidden" in the logarithmic
 *    mean ionization potential.  Much more work is required before this topic
 *    can be fully understood.
 *  - #SSWITCH_LE : Relativistic shell correction.  The Leung, or
 *    relativistic shell correction is a small effect which is due to
 *    relativistic inner shell electrons in very heavy targets.
 *    See Leung, \cite art_ptl1, and Leung, \cite art_ptl2.  #SSWITCH_LE
 *    has no effect unless #SSWITCH_SH is also turned on.
 *  - The Lindhard-Sørensen effect (see lindhard()) is turned on by default.
 *    The Bloch, Mott \& Ahlen effects are included for historical interestest.
 *    Right now these can be turned on by uncommenting a particular section
 *    of the code.
 *  - #SSWITCH_KI : Ultrarelativistic kinematic correction.
 *    This an estimate of the ultrarelativistic kinematic correction from
 *    Ahlen, \cite art_spa2. It corrects to the
 *    finite mass (as opposed to size) of the nucleus in relativistic
 *    electron-nucleus collisions.
 *  - #SSWITCH_RA : Radiative correction.
 *    This is the radiative correction discussed in Ahlen, \cite art_spa2.
 *    It arises from bremsstrahlung
 *    of scattered electrons in ultrarelativistic collisions.  The
 *    form here is that of Jankus, \cite art_vzj.
 *    The parameter Q from that paper is here set equal to the geometric
 *    mean between the the electron rest energy and \f$ 2 m_e c^2 \gamma \f$.
 *  - #SSWITCH_PA : Slowing due to pair production. This value and the value for
 *    the bremsstrahlung correction below are based on the work of
 *    Sørensen, \cite coll_ahs.
 *  - #SSWITCH_BR : Slowing due to projectile bremsstrahlung.  This version is
 *    that of Sørensen, \cite coll_ahs, who has shown that this effect
 *    is much smaller than the version suggested by Weaver \& Westphal,
 *    \cite art_baw3..  This is due to their treatment of the projectile and
 *    target nuclei as a point particles.  That version appeared in some
 *    much older versions of this code, but has been replaced with
 *    Sørensen's version.  We have not yet updated this code to reflect
 *    Sørensen's more recent paper \cite art_ahs1.
 *  - #SSWITCH_BA : Barkas effect.
 *    This is the Barkas correction as calculated in Jackson \&
 *    McCarthy, \cite art_jdj.  It is multiplied
 *    by a factor of two to bring it into agreement with Lindhard, \cite art_jl1.
 *    It is not, however, equal to the
 *    results of Lindhard, and more work is needed to decide which, if any,
 *    form is correct.  The recommended value seems to be the Jackson
 *    \& McCarthy result multiplied by two.  Jackson \& McCarthy do not
 *    have reliable values of \f$ F(V) \f$ for \f$ V < 0.8 \f$ .  For the purposes of the
 *    computation, the cut-off is placed at \f$ V=1.0 \f$ .  I have followed the
 *    convention of Salamon in having the Barkas correction multiply just
 *    the "Bethe" portion of the stopping logarithm rather than the whole
 *    stopping logarithm.  As there is considerable disagreement in the
 *    literature about the application of correction, and as changing
 *    the convention makes makes a difference of less than 1 A MeV even
 *    in calculating the energy of stopping uranium, I have chosen to
 *    leave it where it is.  Furthermore, I have found that a simple
 *    power law \f$ V^{-2} \f$ is adequate to model Jackson \& McCarthy's function
 *    for \f$ V > 1.0 \f$ , so I have used this instead of the numbers found by
 *    reading off one of Jackson \& McCarthy's figures (these values are stored
 *    in the array fva[10], but only the last value is used).
 *
 * \param e1 The projectile kinetic energy in A MeV.
 * \param rel0 Restricted energy loss parameter in eV.
 * \param z0 The projectile charge.
 * \param a1 The projectile atomic number.
 * \param sswitch The switch bit field.
 * \param target A pointer to a ::TDATA structure.
 *
 * \return dE/dx in units of A MeV g<sup>-1</sup> cm<sup>2</sup>
 */
double dedx( double e1, double rel0, double z0, double a1, short sswitch, tdata *target )
{
    static double fva[10]={0.33,0.078,0.03,0.014,0.0084,
        0.0053,0.0035,0.0025,0.0019,0.0014};
    int i;
    double z2,a2;
    double g,b2,b,z1;
    double delt;
    double etam2,cadj;
    double v,fv;
    double S,REL,f1,f2,f3,f4,f6,f7,f8,f9;
    double Sbr=0.0,Bbr;
    double Spa=0.0,dpa,ldpa,l0,Lpa0,Lpa0s,Lpa1,Lpa;
    /* End declarations */
    g=1.0+e1/ATOMICMASSUNIT;
    delt = ( sswitch & SSWITCH_ND ) ? delta(g,target) : olddelta(g,target);
    b2=1.0-1.0/(g*g);
    b=sqrt(b2);
    z2=target->z2;
    a2=target->a2;
    z1 = effective_charge(z0, e1, z2, sswitch);
    f1=0.3070722*z1*z1*z2/(b2*a1*a2);
    f2=log(2.0*ELECTRONMASS*b2/target->iadj);
    if( sswitch & SSWITCH_SH ){
        etam2=1.0/(b*b*g*g);
        cadj=1.0e-6*(target->iadj)*(target->iadj)*etam2*(0.422377
            +etam2*(0.0304043-etam2*0.00038106))
            +1.0e-9*(target->iadj)*(target->iadj)*(target->iadj)*etam2*(3.858019
            +etam2*(-0.1667989 + etam2*0.00157955));
        f2-=cadj/(z2);
        if( sswitch & SSWITCH_LE ){
            f2-=(5.0/3.0)*log(2.0*ELECTRONMASS*b2/target->iadj)*(1.0e+03*target->bind/(z2*ELECTRONMASS))-(target->iadj*target->iadj/(4.0*ELECTRONMASS*ELECTRONMASS*b2));
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
    f8=( sswitch & SSWITCH_KI ) ?
        0.5*(-log(1.0+2.0*((5.4858e-04)*g/a1)) - ((5.4858e-04)*g/a1)*b2/(g*g)) :
        0.0;
    f9=( sswitch & SSWITCH_RA ) ?
        (ALPHA/M_PI)*b2*(6.0822 + log(2.0*g)*(log(2.0*g)*(2.4167 + 0.3333*log(2.0*g))-8.0314)) :
        0.0;
    if( sswitch & SSWITCH_PA ){
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
    if( sswitch & SSWITCH_BR ){
        Bbr=log(1.0 + 2.0*g*0.179524783764566/(exp((1.0/3.0)*log(a1)) + exp((1.0/3.0)*log(a2)))/a1);
        Sbr=5.21721169334564e-07*(z1*z1/a1)*(z1*z1/a1)*(z2*z2/a2)*g*Bbr;
    }
    if( sswitch & SSWITCH_BA ){
        /*
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
    }
    S=f1*(f2*f4+f3+f6-(delt/2.0)+f8+f9) + Sbr + Spa;
    /*
     * Compute restricted energy loss.  REL is activated by setting rel0 >0.
     */
    if(rel0 > 0.0) {
        f7=log(2.0*ELECTRONMASS*b2*g*g/rel0)+b2*(rel0/(2.0*ELECTRONMASS*b2*g*g)-1.0);
        REL=f1*(f2*f4+f3+f6-delt/2.0 - 0.5*f7 +f8);
        return(REL);
    } else {
        return(S);
    }
}
/**
 * \brief Computes the density effect.
 *
 * This function implements the density effect correction as formulated
 * in Sternheimer \& Peierls, \cite art_rms2 and as
 * extended in Sternheimer, Berger \& Seltzer, \cite art_rms1.
 * This version can distinguish between solids and gasses, and between
 * metals and insulators.  For conducting materials, there is a
 * low-energy density effect.
 *
 * \param g Projectile Lorentz factor.
 * \param target A pointer to a ::TDATA structure.
 *
 * \return The value of the density effect.
 */
double delta( double g, tdata *target )
{
    double b,cbar,X,X0,X1;
    /* End declarations */
    X0=target->X0;
    X1=target->X1;
    cbar=2.0*log(target->iadj/target->pla)+1.0;
    b=sqrt(1.0 - 1.0/(g*g));
    X=log10(b*g);
    if(target->etad>0) {
        cbar-=2.303*log10(target->etad);
        X1-=0.5*log10(target->etad);
        X0-=0.5*log10(target->etad);
    }
    if(X < X0) {
        return(target->d0*exp(4.6052*(X-X0)));
    } else if (X >= X0 && X < X1) {
        return(4.6052*X + exp(log(target->a) + target->m*log(X1-X)) - cbar);
    } else {
        return( 4.6052*X - cbar );
    }
}
/**
 * \brief Computes an obsolete version of the density effect.
 *
 * This function implements the density effect correction as originally
 * formulated in Sternheimer \& Peierls, \cite art_rms2.
 * Although it is now obsolete, I have included it here for
 * compatibility with earlier codes.
 *
 * \param g Projectile Lorentz factor.
 * \param target A pointer to a ::TDATA structure.
 *
 * \return The value of the density effect.
 */
double olddelta( double g, tdata *target )
{
    double b,cbar,y;
    double y0,y1,dy3,a;
    /* End declarations */
    if( g < 1.8 ) return(0.0);
    cbar=2.0*log(target->iadj/target->pla)+1.0;
    b=sqrt(1.0 - 1.0/(g*g));
    y=2.0*log(b*g);
    if( target->etad > 0 ){
        y+=log(target->etad);
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
        if(target->iadj>=100.0){
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
/**
 * \brief Computes the Bloch, Mott and Ahlen corrections.
 *
 * This function computes the Mott correction of Ahlen, \cite art_spa1,
 * the Bloch correction of F. Bloch, \cite art_fb1,
 * and the Ahlen correction of Ahlen, \cite art_spa3.
 * All three of these corrections are
 * rendered obsolete by the Lindhard-Sørensen correction, and are
 * included here for historical interest and comparison with older
 * calculations.
 *
 * \param z1 The projectile charge.
 * \param b The projectile velocity in units of the speed of light
 * (\em i.e. \f$ \beta = v/c \f$).
 *
 * \return The sum of the Bloch, Mott and Ahlen corrections.
 *
 * \note The variables lambda and theta0 are
 * free parameters in the Ahlen correction.  Theta0 also appears in
 * the Mott correction.  Here I have used Ahlen's recommended values,
 * lambda = 1, theta0 = 0.1.
 * An alternative formula, \f$ \theta_0 = \sqrt{\alpha/(\beta \gamma \lambda)} \f$ , is
 * suggested by Waddington, Freier \& Fixsen, \cite art_cjw1.
 *
 * \warning The Mott correction has a severely
 * limited range of validity, especially for high charges.  It's so
 * bad it can render the calculation not just inaccurate, but
 * unphysical (dE/dx \< 0) below about 10 A MeV for uranium.  Ahlen
 * recommends turning the Mott correction off for \f$ Z/\beta > 100 \f$.
 * Here for \f$ Z/\beta > 100 \f$ the Mott correction is given the value at
 * \f$ Z/\beta = 100 \f$. This prescription is given by Waddington,
 * Freier \& Fixsen, \cite art_cjw1.
 *
 * \bug Currently, this function is not called by anything.
 */
double bma( double z1, double b )
{
    int n,msum;
    double b2,y,y2,sumr,fn,fn2,f3,f5;
    double lambda, theta0, cosx, st;
    gsl_complex Cz1,Cz2;
    /* End declarations */
    /*
     * Compute a sum needed by the Bloch Correction
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
     * Compute the Bloch and Ahlen Corrections
     */
    lambda=1.0;
    theta0=0.1;
    f3=-y2*(1.202+sumr) + relbloch(z1,b,lambda,theta0);
    Cz1=gsl_complex_rect(0.5,-y);
    Cz2=gsl_complex_rect(1.0,y);
    cosx=cos(2.0*(GSL_IMAG(complex_lngamma(Cz1))+GSL_IMAG(complex_lngamma(Cz2))));
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
/**
 * \brief Compute the relativistic Bloch correction.
 *
 * This is the relativistic Bloch (or Ahlen) correction of Ahlen,
 * \cite art_spa3.  The evaluation of this correction has
 * been enormously simplified by the use of fully complex arithmetic.
 *
 * \param z12 The projectile charge.
 * \param b1 The projectile velocity in units of the speed of light
 * (\em i.e. \f$ \beta = v/c \f$).
 * \param lambda A free parameter, described in bma().
 * \param theta0 A free parameter, described in bma().
 *
 * \return The value of the relativistic Bloch correction.
 */
double relbloch( double z12, double b1, double lambda, double theta0 )
{
    gsl_complex Cl1, Cl2;
    gsl_complex Cf1,Cf2,Cf3,Cf4,Cf5,Cf6,Cf7;
    gsl_complex Ci,Cnu,Cth,Cabgl;
    double ge = 0.5772157;
    double abgl;
    double g,nu;
    double sigma,es,bloch2;
    /* End declarations */
    nu=z12*ALPHA/b1;
    g=1.0/sqrt(1.0-b1*b1);
    abgl=ALPHA/(b1*g*lambda);
    sigma=GSL_IMAG(complex_lngamma(gsl_complex_rect(1.0,nu)));
    Ci=gsl_complex_rect(0.0,1.0);
    Cnu=gsl_complex_rect(1.0, 2.0*nu);
    Cth=gsl_complex_rect(theta0/2.0,0.0);
    Cabgl=gsl_complex_rect(abgl,0.0);
    Cf1=gsl_complex_div(gsl_complex_rect(0.0,-1.0),Cnu);
    Cf2=gsl_complex_pow(Cth,Cnu);
    Cf3=gsl_complex_pow(Cabgl,Cnu);
    Cf4=gsl_complex_exp(gsl_complex_rect(0.0,2.0*sigma));
    Cf5=gsl_complex_div(gsl_complex_rect(1.0,0.0),Cnu);
    Cf6=gsl_complex_sub(gsl_complex_rect(log(theta0/2.0),0.0),Cf5);
    Cf7=gsl_complex_add(gsl_complex_rect(log(4.0/abgl)+ge-1.0,0.0),Cf5);
    Cl1=gsl_complex_mul(Cf1,
        gsl_complex_sub(gsl_complex_mul_real(Cf2,2.0),gsl_complex_mul(Cf3,Cf4)));
    Cl2=gsl_complex_mul(Cf1,
        gsl_complex_add(
            gsl_complex_mul_real(
                gsl_complex_mul(Cf2,Cf6),2.0),gsl_complex_mul(Cf3,
            gsl_complex_mul(Cf4,Cf7))));
    es=(M_PI*nu)*exp(M_PI*nu)/sinh(M_PI*nu);
    bloch2=(M_PI/2.0)*b1*b1*nu*es*(4.0*nu*log(2.0)*GSL_REAL(Cl1)
        + (nu*M_PI-1.0)*GSL_IMAG(Cl1)
        + 2.0*nu*GSL_REAL(Cl2));
    return(bloch2);
}
/**
 * \brief Compute the Lindhard-Sørensen correction.
 *
 * This is the Lindhard-Sørensen correction including finite nuclear
 * size effects as described in Lindhard \& Sørensen, \cite art_jl2.
 * The defined variable #SSWITCH_NS will turn off the
 * nuclear size effect if it is set to zero.  For values of the Lorentz
 * factor above 10/R, where R is the nuclear size divided by the electron
 * Compton wavelength, the correction is set to its asymptotic value
 * which is described by Sørensen, \cite proc_ahs. This also avoids some
 * difficulties with the evaluation of the confluent hypergeometric function
 * (A. H. Sørensen, private communication).
 *
 * \param zz The projectile charge.
 * \param aa The projectile atomic mass.
 * \param bb The projectile velocity in units of the speed of light
 * (\em i.e. \f$ \beta = v/c \f$).
 * \param sswitch The switch bit field.
 *
 * \return The value of the Lindhard-Sørensen correction.
 */
double lindhard( double zz, double aa, double bb, short sswitch )
{
    gsl_complex Cexir, Cexis, Cedr, Ceds, Cske, Cmske, Cgrgs;
    gsl_complex Clamr,Clams;
    gsl_complex Caar, Cbbr, Caas, Cbbs, Czzr;
    gsl_complex Cpi,Cone;
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
    /* End declarations */
    Cpi=gsl_complex_rect(M_PI,0.0);
    Cone=GSL_COMPLEX_ONE;
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
    if ((gg < 10.0/rho) || !(sswitch & SSWITCH_NS)) {
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
                Cske=gsl_complex_rect(sk+1.0,eta);
                Cexir=gsl_complex_sqrt(gsl_complex_div(
                    gsl_complex_rect(k,-eta/gg),gsl_complex_rect(sk,-eta)));
                Cedr=gsl_complex_mul(Cexir,
                    gsl_complex_exp(gsl_complex_rect(0.0,
                        -GSL_IMAG(complex_lngamma(Cske))+(M_PI/2.0)*(l-sk))));
                if( sswitch & SSWITCH_NS ){
                    Cmske=gsl_complex_rect(-sk+1.0,eta);
                    Cexis=gsl_complex_sqrt(gsl_complex_div(gsl_complex_rect(k,-eta/gg),gsl_complex_rect(-sk,-eta)));
                    Ceds=gsl_complex_mul(Cexis,
                        gsl_complex_exp(gsl_complex_rect(0.0,
                            -GSL_IMAG(complex_lngamma(Cmske)) + (M_PI/2.0)*(l+sk))));
                    Caar=Cske;
                    Caas=Cmske;
                    Cbbr=gsl_complex_rect(2.0*sk + 1.0,0.0);
                    Cbbs=gsl_complex_rect(-2.0*sk + 1.0,0.0);
                    Czzr=gsl_complex_rect(0.0,2.0*prh);
                    Clamr=gsl_complex_mul(gsl_complex_mul(Cexir,gsl_complex_exp(gsl_complex_rect(0.0,-prh))),
                        complex_hyperg(Caar,Cbbr,Czzr));
                    Clams=gsl_complex_mul(gsl_complex_mul(Cexis,gsl_complex_exp(gsl_complex_rect(0.0,-prh))),
                        complex_hyperg(Caas,Cbbs,Czzr));
                    grgs=GSL_IMAG(Clamr)/GSL_IMAG(Clams);
                    Cgrgs=complex_lngamma(Cbbs);
                    grgs*=exp( GSL_REAL(complex_lngamma(Cske)) - GSL_REAL(complex_lngamma(Cmske))
                        -GSL_REAL(complex_lngamma(Cbbr)) + GSL_REAL(Cgrgs)
                        +2.0*sk*log(2.0*prh));
                    if(cos(GSL_IMAG(Cgrgs))<1.0)grgs*= -1.0;
                    if(fabs(grgs) > 1.0e-9) {
                        frgr=sqrt((gg-1.0)/(gg+1.0))*GSL_REAL(Clamr)/GSL_IMAG(Clamr);
                        fsgs=sqrt((gg-1.0)/(gg+1.0))*GSL_REAL(Clams)/GSL_IMAG(Clams);
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
                    Ceds=GSL_COMPLEX_ZERO;
                }
                r=1.0 + H*H + 2.0*H*(GSL_REAL(Cedr)*GSL_REAL(Ceds) + GSL_IMAG(Cedr)*GSL_IMAG(Cedr));
                dk[i]=gsl_complex_arg(gsl_complex_add(Cedr,gsl_complex_mul_real(Ceds,H)));
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
/**
 * \brief Compute a mathematical function related to bremsstrahlung.
 *
 * This function is used in an obsolete version of projectile slowing
 * due to nuclear-nuclear bremsstrahlung.  It appears in Heitler's treatment
 * of bremsstrahlung, \cite book_wh, which was adapted by
 * Weaver \& Westphal, \cite art_baw3.
 *
 * \param x The input parameter.
 *
 * \return The value of the function.
 *
 * \bug Currently, this function is unused.
 */
double Fbrems( double x )
{
    int n=1;
    double t=1.0,s=0.0;
    /* End declarations */
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
/**
 * \brief Computes total range given initial energy.
 *
 * This function computes total range given initial energy.  The technique
 * is quite clever, in that if from one call to the next, the projectile
 * and target material parameters do not change, the calculation of
 * range is performed by table interpolation rather than direct integration.
 * The savings in calculation time can be enormous.  However, the
 * range of valid energies is limited by the size of the table.  The
 * function dE/dx is evaluated at most of the energies defined by the function
 * energy_table().  Results are stored in the ::RANGE_TABLE array
 * ::trange.  The size of ::trange (set by ::MAXAB) is arbitrary
 * and should be set to whatever is most useful.  Certainly it should
 * be no smaller than the number of target materials being used.
 *
 * \param e Initial projectile kinetic energy in A MeV.
 * \param z1 Projectile charge.
 * \param a1 Projectile atomic mass.
 * \param sswitch The switch bit field.
 * \param target A pointer to a ::TDATA structure.
 * \param tno A pointer to the index of the most recently generated or used table.
 *
 * \return Projectile range in g cm<sup>-2</sup>.
 */
double range( double e, double z1, double a1, short sswitch, tdata *target, int *tno )
{
    extern range_table trange[MAXAB];
    int table=1,i,k,l;
    double de2,dr,dedx1,dedx2,dedx3,dedx4,e0,e1,e2,e3,e4,dt,dtk;
    double rel = 0.0;
    time_t ct;
    /* End declarations */
    /*
     * Search the range table for existing data
     */
    for (k=0;k<MAXAB;k++) {
        if ((trange[k].z1 == z1) && (trange[k].a1 == a1) &&
            (trange[k].sswitch == sswitch) &&
            strncmp(trange[k].target->name,target->name,NAMEWIDTH) == 0) {
            *tno = k;
            table = 0;
            break;
        }
    }
    if(table){
        /*
         * Define the meta-data for the new table.  First, find an unused index
         * in the trange array.
         */
        k=0;
        while (k < MAXAB && trange[k].target != NULL) k++;
        /*
         * If all the tables are filled, replace the oldest entry.
         */
        if (k == MAXAB) {
            l=0;
            ct = time(NULL);
            dt = 0.0;
            while (l < MAXAB) {
                dtk = difftime(ct,trange[l].timestamp);
                if (dtk > dt) {
                    dt = dtk;
                    k = l;
                }
                l++;
            }
        }
        *tno=k;
        trange[*tno].z1 = z1;
        trange[*tno].a1 = a1;
        trange[*tno].sswitch = sswitch;
        trange[*tno].target = target;
        trange[*tno].timestamp = time(NULL);
        i=0;
        while(energy_table(i) < 8.0) {
            e1=energy_table(i);
            trange[*tno].range[i]=benton(e1, z1, a1, target);
            i++;
        }
        while(i<MAXE){
            e0=energy_table(i-1);
            de2=(energy_table(i)-e0)/2.0;
            e1=e0+1.33998104*de2;
            dedx1=dedx(e1,rel,z1,a1,sswitch,target);
            e2=e0+1.86113631*de2;
            dedx2=dedx(e2,rel,z1,a1,sswitch,target);
            e3=e0+0.13886369*de2;
            dedx3=dedx(e3,rel,z1,a1,sswitch,target);
            e4=e0+0.66001869*de2;
            dedx4=dedx(e4,rel,z1,a1,sswitch,target);
            dr=de2*(0.65214515/dedx1 + 0.34785485/dedx2 + 0.34785485/dedx3
                + 0.65214515/dedx4);
            trange[*tno].range[i] = trange[*tno].range[i-1] + dr;
            i++;
        }
    }
    if( e > energy_table(0)) {
        i=1;
        while( e > energy_table(i) ) i++;
        return( trange[*tno].range[i-1] + ( e - energy_table(i-1) )
            *(trange[*tno].range[i]-trange[*tno].range[i-1])/(energy_table(i)-energy_table(i-1)) );
    } else {
        return( e*trange[*tno].range[0]/energy_table(0) );
    }
}
/**
 * \brief Computes total range by direct integration of dE/dx.
 *
 * This function computes total range by direct integration of the
 * dedx() function.  It does not create a range table or do table
 * interpolation.
 *
 * \param e Initial energy in A MeV.
 * \param z1 Projectile charge.
 * \param a1 Projectile mass.
 * \param sswitch The switch bit field.
 * \param target A pointer to a ::TDATA structure.
 *
 * \return Total range in g cm<sup>-2</sup>.
 * \bug Currently, this function isn't called by anything.
 */
double qrange( double e, double z1, double a1, short sswitch, tdata *target )
{
    int i;
    double de2,dr,dedx1,dedx2,dedx3,dedx4,e1,e2,e3,e4;
    double ei=8.0,lef,lei,en=1.0,pen,rel = 0.0;
    double entries;
    double ra[MAXE];
    /* End declarations */
    if(e > ei){
        ra[0]=benton(ei,z1,a1,target);
        i=1;
        lei=log10(ei);
        lef=log10(e);
        entries = (double) (MAXE - 1);
        do {
            pen=en;
            en=exp(M_LN10*(lei + ((double)i)*lef/entries));
            de2=(en-pen)/2.0;
            e1=pen+1.33998104*de2;
            dedx1=dedx(e1,rel,z1,a1,sswitch,target);
            e2=pen+1.86113631*de2;
            dedx2=dedx(e2,rel,z1,a1,sswitch,target);
            e3=pen+0.13886369*de2;
            dedx3=dedx(e3,rel,z1,a1,sswitch,target);
            e4=pen+0.66001869*de2;
            dedx4=dedx(e4,rel,z1,a1,sswitch,target);
            dr=de2*(0.65214515/dedx1 + 0.34785485/dedx2 + 0.34785485/dedx3
                + 0.65214515/dedx4);
            ra[i] = ra[i-1] + dr;
            i++;
        }while(en < e);
        return(ra[i-1]);
    } else if (e > 1.0 && e <= ei) {
        return( benton(e,z1,a1,target) );
    } else {
        return( e*benton(en,z1,a1,target)/en );
    }
}
/**
 * \brief Computes ranges at low energies.
 *
 * This function is the result of empirical fits to very low energy
 * 1 A MeV \< E \< 8 A MeV ion ranges.  It follows the methods of
 * Barkas \& Berger, \cite coll_whb. A simplified
 * discussion, with a more complicated formula is given in
 * Benton \& Henke, \cite art_evb1.
 * As yet I know of no nicer way to deal with these low energies.
 *
 * \param e Projectile kinetic energy in A MeV.
 * \param z1 Projectile charge.
 * \param a1 Projectile atomic mass.
 * \param target A pointer to a ::TDATA structure.
 *
 * \return Projectile range in g cm<sup>-2</sup>.
 *
 * \note The array join[4] demarcates three energy regions represented by
 * the three sets of coefficients in amn[3][4][4]. The demarcation is variable
 * in order to minimize discontinuities at the boundary.  The coefficients
 * in cjoin[2][7], which is used to initialize join[4], are inherited from
 * legacy code; I have not found them in the non-obscure literature.
 * Approximately, the three regions are \em E \< 1 A MeV, 1 \< \em E \< 7 A MeV
 * and \em E \> 7 A MeV.  I can find no reason why join[4] has four elements and
 * not two.
 */
double benton( double e, double z1, double a1, tdata *target )
{
    int l,m,n;
    double g,b,bzz,x;
    double term,logt,logi,loglambda;
    double cr = PROTONMASS/ATOMICMASSUNIT;
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
    /* End declarations */
    logt=log( e * cr );
    logi=log( target->iadj );
    for (l = 0; l < 2; l++) {
        join[l] = cjoin[l][m=6];
        while (m > 0) join[l] = 0.001*target->iadj*join[l] + cjoin[l][--m];
    }
    for (l = 0; l < 3; l++) {
        loglambda=0.0;
        for (m = 3; m >= 0; m-- ) {
            term=amn[l][m][n=3];
            while (n > 0) term=term*logt + amn[l][m][--n];
            loglambda+=term;
            loglambda*=logi;
        }
        loglambda/=logi;
        loglambda+=log( (target->a2)/(target->z2) );
        prnglo[l]=exp(loglambda);
        if (l == 1) prnglo[l]*=1.0e-03;
    }
    g=1.0 + e/ATOMICMASSUNIT;
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
        *( (target->a2)/(target->z2) )*1.0e-06*exp((8.0/3.0)*log( z1 ));
    return(( (a1/cr)/(z1*z1) )*(prnglo[l] + bzz*cz[n]));
}
/**
 * \brief Extract energies from range tables.
 *
 * This function extracts energies from a range table by table interpolation.
 * It calls range() to initialze the range table or to find the correct
 * table if it has already been computed.
 *
 * \param e Projectile kinetic energy [A MeV].
 * \param r0 Range [g cm<sup>-2</sup>].
 * \param z1 Projectile charge.
 * \param a1 Projectile atomic mass.
 * \param sswitch The switch bit field.
 * \param target A pointer to a ::TDATA structure.
 *
 * \return The final energy of the projectile.
 */
double renergy( double e, double r0, double z1, double a1, short sswitch, tdata *target )
{
    extern range_table trange[MAXAB];
    double rr,r;
    int i,k;
    int tno = 0;
    /* End declarations */
    if( e > 0.0 ){
        rr = range(e,z1,a1,sswitch,target,&tno);
        r = rr - r0;
        if( r < 0.0 ) return(0.0);
    } else {
        rr = range(energy_table(0),z1,a1,sswitch,target,&tno);
        r = r0;
    }
    if(r > trange[tno].range[0]){
        i=1;
        while( r > trange[tno].range[i] ) i++;
        return(energy_table(i-1)+(r-trange[tno].range[i-1])*(energy_table(i)-energy_table(i-1))/
            (trange[tno].range[i]-trange[tno].range[i-1]));
    } else {
        return(energy_table(0)*r/trange[tno].range[0]);
    }
}
/**
 * \brief Parses and executes the task list.
 *
 * This utility function steps through the range, energy and dE/dx tasks
 * specified in the input data file.
 * The tasks are denoted by a single letter:
 *     - \em r compute ranges
 *     - \em e compute energies
 *     - \em d compute dE/dx
 *     - \em j compute dJ/dx (primary ionization)
 *
 * The task letter should be followed by the energy (or range) at which to
 * compute range (or energy), the charge and mass of the particle, and the
 * name of the target material.  Names of target materials can be found in the
 * target.ini file.  Target material names may be up to #NAMEWIDTH
 * characters in length and should contain no whitespace.
 *
 * \param finput An open file pointer containing the task list.
 * \param foutput An open file pointer to write results to.
 * \param sswitch The switch bit field.
 * \param extratargets A pointer to an array of ::TDATA structures.
 *
 * \bug The primary ionization parameters are currently hard-coded.
 */
void run_range( FILE *finput, FILE *foutput, short sswitch, tdata *extratargets )
{
    tdata *target;
    char task[2];
    char abs[NAMEWIDTH+1];
    double red1,red2,z1,a1;
    double out=0.0;
    int icols=6;
    int k=0;
    int tno=0;
    /* End declarations */
    for(;;){
        icols=fscanf(finput,"%s %lf %lf %lf %lf %s\n",
            task,&red1,&red2,&z1,&a1,abs);
        if(icols < 6) break;
        target = find_target(abs,extratargets);
        out=0.0;
        if(icols==6 && strncmp( target->name, "Unknown", NAMEWIDTH ) != 0){
            if(strncmp( task, "r", 1 )==0){
                out=range(red1,z1,a1,sswitch,target,&tno);
            } else if(strncmp( task, "e", 1 )==0){
                out=renergy(red1,red2,z1,a1,sswitch,target);
            } else if(strncmp( task, "d", 1 )==0){
                out=dedx(red1,red2,z1,a1,sswitch,target);
            } else if(strncmp( task, "j", 1 )==0){
                out=djdx(red1,z1,2.0,0.05,3.04,sswitch,target);
            }
        }
        fprintf(foutput,"%f\n",out);
    }
    return;
}
/**
 * \brief Initializes the value of of the switch bit field.
 *
 * This utility reads an INI-type file and sets the switch bit field
 * accordingly.
 *
 * \param switchfile The name of an INI-type file containing switch configuration.
 *
 * \return The switch bit field.
 *
 * \warning If the iniparser library is not found, this function will only
 * return the default value #SSWITCH_DEFAULT.
 */
short init_switch( char *switchfile )
{
    short sswitch = SSWITCH_DEFAULT;
#ifdef HAVE_INIPARSER_H
    dictionary *d;
    /* End declarations */
    if (access(switchfile,R_OK) == 0) {
        sswitch = 0;
        d = iniparser_load( switchfile );
        if (iniparser_getboolean(d,"switch:barkas",0)) sswitch |= SSWITCH_BA;
        if (iniparser_getboolean(d,"switch:shell",0))  sswitch |= SSWITCH_SH;
        if (iniparser_getboolean(d,"switch:leung",0))  sswitch |= SSWITCH_LE;
        if (iniparser_getboolean(d,"switch:new delta",1))  sswitch |= SSWITCH_ND; /* True by default! */
        if (iniparser_getboolean(d,"switch:new electron capture",0))  sswitch |= SSWITCH_EC;
        if (iniparser_getboolean(d,"switch:finite nuclear size",1))  sswitch |= SSWITCH_NS; /* True by default! */
        if (iniparser_getboolean(d,"switch:kinematic",0))  sswitch |= SSWITCH_KI;
        if (iniparser_getboolean(d,"switch:radiative",0))  sswitch |= SSWITCH_RA;
        if (iniparser_getboolean(d,"switch:pair",0))  sswitch |= SSWITCH_PA;
        if (iniparser_getboolean(d,"switch:bremsstrahlung",0))  sswitch |= SSWITCH_BR;
        iniparser_freedict(d);
    } else {
        fprintf(stderr,"Could not read switch file: %s. Using default crange switches.\n",switchfile);
    }
#endif
    return(sswitch);
}
/**
 * \brief Read optional target data file.
 *
 * This utility reads an INI-type file and returns an array of pointers to
 * ::TDATA structures.
 *
 * \param targetfile the name of an INI-type file containing target data.
 *
 * \return A pointer to an array of ::TDATA structures.  This pointer must
 * be free()d!
 *
 * \warning If the iniparser library is not found, this function will only
 * return a NULL pointer.
 */
tdata *init_target( char *targetfile )
{
    tdata *table;
    table = NULL;
#ifdef HAVE_INIPARSER_H
    dictionary *d;
    char *section;
    char key[NAMEWIDTH+5+1];
    int i,n;
    /* End declarations */
    if (access(targetfile,R_OK) == 0) {
        d = iniparser_load( targetfile );
        n = iniparser_getnsec(d);
        table = (tdata*)calloc(n+1,sizeof(tdata));
        for (i=0;i<n;i++) {
            section=iniparser_getsecname(d,i);
            sprintf(key,"%s:%s",section,"name");
            /* printf("%s = %s\n",key,iniparser_getstring(d,key,"Unknown")); */
            strncpy(table[i].name,iniparser_getstring(d,key,"Unknown"),NAMEWIDTH);
            sprintf(key,"%s:%s",section,"z2");
            /* printf("%s = %f\n",key,iniparser_getdouble(d,key,0.0)); */
            table[i].z2 = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"a2");
            table[i].a2 = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"iadj");
            table[i].iadj = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"rho");
            table[i].rho = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"pla");
            table[i].pla = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"etad");
            table[i].etad = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"bind");
            table[i].bind = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"x0");
            table[i].X0 = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"x1");
            table[i].X1 = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"a");
            table[i].a = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"m");
            table[i].m = iniparser_getdouble(d,key,0.0);
            sprintf(key,"%s:%s",section,"d0");
            table[i].d0 = iniparser_getdouble(d,key,0.0);
        }
        iniparser_freedict(d);
        /*
         * Terminate the table with a dummy value.
         */
        strncpy(table[n].name,"Unknown",NAMEWIDTH);
    } else {
        fprintf(stderr,"Could not read target file: %s. Using built-in crange targets.\n",targetfile);
    }
#endif
    return(table);
}
/**
 * \brief Initialize range-energy tables.
 *
 * Sets all data in the array of ::RANGE_TABLE structures trange to zero.
 *
 */
void init_table(void)
{
    extern range_table trange[MAXAB];
    int k,l;
    /* End declarations */
    for (k=0;k<MAXAB;k++) {
        trange[k].z1 = trange[k].a1 = 0.0;
        trange[k].sswitch = 0;
        trange[k].target = NULL;
        trange[k].timestamp = (time_t)0;
        for (l=0;l<MAXE;l++) trange[k].range[l] = 0.0;
    }
    return;
}
/**
 * \brief Returns the energy corresponding to a value in a range table.
 *
 * This utility returns an energy value from a (virtual) vector containing
 * A logarithmically uniform distribution of energies between a minimum
 * and maximum energy (defined by #LOGTENEMIN and #LOGTENEMAX),
 * with a number of entries given by #MAXE.
 *
 * \param i The index of the vector.
 *
 * \return The \em i th energy in A MeV.
 */
double energy_table( int i )
{
    /*
     * The values of all of these are known at compile time.
     */
    double decades, entries;
    /* End declarations */
    decades=LOGTENEMAX-LOGTENEMIN;
    entries=MAXE - 1.0;
    return(exp(M_LN10*(LOGTENEMIN + ((double)i)*decades/entries)));
}
/**
 * \brief Finds target data corresponding to a target name.
 *
 * This function returns a pointer to a structure containing the target
 * data corresponding to the input name.  There is a built-in list. The
 * built-in list may be added to or overridden by supplying an INI-type
 * file on the command line, which will then be parsed & passed to this
 * function.  If the special target name "List" is passed to this function,
 * the built-in list will be printed as an INI-type file.
 *
 * \param target The name of a target.
 * \param extratargets A pointer to an array of ::TDATA structures.
 *
 * \return A pointer to a structure containing the target data.  If the
 * name of the target was not found, it will point to a dummy structure.
 */
tdata *find_target( char *target, tdata *extratargets )
{
    int k=0;
    static tdata targets[] = {
        /* name     ,     z2,      a2,  iadj,        rho,        pla,etad,       bind,      X0,     X1,       a,      m,   d0*/
        { "H"       ,  1.000,   1.008,  19.2, 8.3748e-05, 2.6300e-01, 1.0, 1.3606e-02,  1.8639, 3.2718, 0.14092, 5.7273, 0.00 },
        { "He"      ,  2.000,   4.003,  41.8, 1.6632e-04, 2.6300e-01, 1.0, 7.7872e-02,  2.2017, 3.6122, 0.13443, 5.8347, 0.00 },
        { "C"       ,  6.000,  12.011,  78.0, 2.0000e+00, 2.8803e+01, 0.0, 1.0251e+00, -0.0351, 2.4860, 0.20240, 3.0036, 0.10 },
        { "N"       ,  7.000,  14.007,  82.0, 1.1653e-03, 6.9500e-01, 1.0, 1.4782e+00,  1.7378, 4.1323, 0.15349, 3.2125, 0.00 },
        { "O"       ,  8.000,  15.999,  95.0, 1.3315e-03, 7.4400e-01, 0.0, 2.0359e+00,  1.7541, 4.3213, 0.11778, 3.2913, 0.00 },
        { "Na"      , 11.000,  22.990, 149.0, 9.7100e-01, 1.9641e+01, 0.0, 4.4097e+00,  0.2880, 3.1962, 0.07772, 3.6452, 0.08 },
        { "Al"      , 13.000,  26.980, 166.0, 2.6989e+00, 3.2860e+01, 0.0, 6.5929e+00,  0.1708, 3.0127, 0.08024, 3.6345, 0.12 },
        { "Si"      , 14.000,  28.090, 173.0, 2.3300e+00, 3.1055e+01, 0.0, 7.8751e+00,  0.2014, 2.8715, 0.14921, 3.2546, 0.14 },
        { "P"       , 15.000,  30.974, 173.0, 2.2000e+00, 2.9743e+01, 0.0, 9.2905e+00,  0.1696, 2.7815, 0.23610, 2.9158, 0.14 },
        { "Ar"      , 18.000,  39.948, 188.0, 1.6620e-03, 7.8900e-01, 1.0, 1.4382e+01,  1.7635, 4.4855, 0.19714, 2.9618, 0.00 },
        { "Fe"      , 26.000,  55.847, 286.0, 7.8740e+00, 5.5172e+01, 0.0, 3.4582e+01, -0.0012, 3.1531, 0.14680, 2.9632, 0.12 },
        { "Ni"      , 28.000,  58.690, 311.0, 8.9020e+00, 5.9385e+01, 0.0, 4.1324e+01, -0.0566, 3.1851, 0.16496, 2.8430, 0.10 },
        { "Cu"      , 29.000,  63.550, 323.0, 8.9600e+00, 5.8270e+01, 0.0, 4.4973e+01, -0.0254, 3.2792, 0.14339, 2.9044, 0.08 },
        { "Ge"      , 32.000,  72.610, 350.0, 5.3230e+00, 4.4141e+01, 0.0, 5.7047e+01,  0.3376, 3.6096, 0.07188, 3.3306, 0.14 },
        { "Ag"      , 47.000, 107.870, 470.0, 1.0500e+01, 6.1635e+01, 0.0, 1.4451e+02,  0.0657, 3.1074, 0.24585, 2.6899, 0.14 },
        { "Ba"      , 56.000, 137.330, 491.0, 3.5000e+00, 3.4425e+01, 0.0, 2.2118e+02,  0.4190, 3.4547, 0.18268, 2.8906, 0.14 },
        { "Os"      , 76.000, 190.200, 746.0, 2.2570e+01, 8.6537e+01, 0.0, 4.6939e+02,  0.0891, 3.5414, 0.12751, 2.9608, 0.10 },
        { "Pt"      , 78.000, 195.090, 790.0, 2.1450e+01, 8.4389e+01, 0.0, 5.0101e+02,  0.1484, 3.6212, 0.11128, 3.0417, 0.12 },
        { "Au"      , 79.000, 196.970, 790.0, 1.9320e+01, 8.0215e+01, 0.0, 5.1732e+02,  0.2021, 3.6979, 0.09756, 3.1101, 0.14 },
        { "Pb"      , 82.000, 207.200, 823.0, 1.1350e+01, 6.1072e+01, 0.0, 5.6834e+02,  0.3776, 3.8073, 0.09359, 3.1608, 0.14 },
        { "U"       , 92.000, 238.030, 890.0, 1.5370e+01, 7.7986e+01, 0.0, 7.6220e+02,  0.2260, 3.3721, 0.19677, 2.8171, 0.14 },
        { "Air"     ,  7.312,  14.667,  85.4, 1.2048e-03, 7.0700e-01, 1.0, 1.7150e+00,  1.7418, 4.2759, 0.10914, 3.3994, 0.00 },
        { "ArCO2"   , 15.867,  34.892, 174.7, 1.7000e-03, 8.0120e-01, 1.0, 1.7150e+00,  1.7418, 4.2759, 0.10914, 3.3994, 0.00 },
        { "BC-408"  ,  3.381,   6.248,  62.8, 1.0320e+00, 2.1534e+01, 1.0, 4.9530e-01,  0.1769, 2.6747, 0.11442, 3.3762, 0.00 },
        { "BP-1"    , 14.795,  32.575, 242.5, 3.0000e+00, 3.3636e+01, 0.0, 2.7489e+01,  0.0843, 3.6297, 0.06445, 3.3655, 0.00 },
        { "CH2"     ,  2.667,   4.676,  57.4, 9.4000e-01, 2.1099e+01, 0.0, 3.5078e-01,  0.1370, 2.5177, 0.12108, 3.4292, 0.00 },
        { "CO2"     ,  7.333,  14.670,  85.0, 1.8421e-03, 8.7400e-01, 1.0, 1.6990e+00,  1.6294, 4.1825, 0.11768, 3.3227, 0.00 },
        { "CR-39"   ,  3.946,   7.413,  73.1, 1.4000e+00, 2.4780e+01, 0.0, 7.2426e-01,  0.1562, 2.6507, 0.12679, 3.3076, 0.00 },
        { "CsI"     , 54.000, 129.905, 553.1, 4.5100e+00, 3.9455e+01, 0.0, 2.0258e+02,  0.0395, 3.3353, 0.25381, 2.6657, 0.00 },
        { "Halo"    ,  1.077,   1.238,  19.2, 4.0880e-24, 5.4342e-11, 1.0, 1.8554e-02,  2.0000, 3.5000, 0.13500, 5.7500, 0.00 },
        { "Hosta"   ,  4.250,   8.091,  72.3, 1.4000e+00, 2.4722e-01, 0.0, 7.7212e-01,  0.1606, 2.6255, 0.12860, 3.3288, 0.00 },
        { "ISM"     ,  1.077,   1.238,  19.2, 2.0440e-24, 3.8426e-11, 1.0, 1.8554e-02,  2.0000, 3.5000, 0.13500, 5.7500, 0.00 },
        { "Kapton"  ,  5.026,   9.803,  75.9, 1.4200e+00, 2.4586e+01, 0.0, 9.1859e-01,  0.1509, 2.5631, 0.15972, 3.1912, 0.00 },
        { "Kevlar"  ,  4.000,   7.567,  71.7, 1.4500e+00, 2.5229e+01, 0.0, 9.1859e-01,  0.1509, 2.5631, 0.15972, 3.1912, 0.00 },
        { "Lexan"   ,  4.061,   7.706,  73.1, 1.2040e+00, 2.2915e+01, 0.0, 6.8789e-01,  0.1606, 2.6255, 0.12860, 3.3288, 0.00 },
        { "LH2"     ,  1.000,   1.008,  21.8, 6.0000e-02, 7.0310e+00, 0.0, 1.3606e-02,  0.4759, 1.9215, 0.13483, 5.6249, 0.00 },
        { "Mesh"    , 21.400,  44.200, 223.0, 2.1500e+00, 2.9400e+01, 0.0, 3.4582e+01, -0.0012, 3.1531, 0.14680, 2.9632, 0.12 },
        { "Mylar"   ,  4.456,   8.735,  78.7, 1.4000e+00, 2.4595e+01, 0.0, 8.4108e-01,  0.1562, 2.6507, 0.12679, 3.3067, 0.00 },
        { "SiO2"    , 10.000,  20.029, 139.2, 2.3200e+00, 3.1014e+01, 0.0, 3.9822e+00,  0.1385, 3.0025, 0.08408, 3.5064, 0.00 },
        { "Teflon"  ,  8.000,  16.669,  99.1, 2.2000e+00, 2.9609e+01, 0.0, 2.1465e+00,  0.1648, 2.7404, 0.10606, 3.4046, 0.00 },
        { "Water"   ,  3.333,   6.005,  75.0, 1.0000e+00, 2.1469e+01, 0.0, 6.8770e-01,  0.2400, 2.8004, 0.09116, 3.4773, 0.00 },
        /* THIS MUST BE THE LAST STRUCTURE DEFINITION! */
        { "Unknown" ,  0.000,   0.000,   0.0, 0.0000e+00, 0.0000e+00, 0.0, 0.0000e+00,  0.0000, 0.0000, 0.00000, 0.0000, 0.00 }
    };
    /* End declarations */
    if (strncmp("List", target, NAMEWIDTH) == 0) {
        while (strncmp(targets[k].name, "Unknown", NAMEWIDTH) != 0) {
            print_target(&targets[k]);
            k++;
        }
        print_target(&targets[k]);
        return(&targets[0]);
    }
    if (extratargets != NULL) {
        while (strncmp(extratargets[k].name, "Unknown", NAMEWIDTH) != 0) {
            if (strncmp(extratargets[k].name, target, NAMEWIDTH) == 0) return(&extratargets[k]);
            k++;
        }
    }
    k=0;
    while (strncmp(targets[k].name, "Unknown", NAMEWIDTH) != 0) {
        if (strncmp(targets[k].name, target, NAMEWIDTH) == 0) break;
        k++;
    }
    /*
     * If the while loop doesn't break, the 'Unknown' target will be returned.
     */
    return(&targets[k]);
}
/**
 * \brief Prints a target table entry in INI format.
 *
 * This utility prints a ::TDATA structure in INI format.
 *
 * \param target A pointer to a ::TDATA structure.
 */
void print_target( tdata *target )
{
    printf("[%s]\n",target->name);
    printf("name = %10s ; Target name\n", target->name);
    printf("z2   = %10.3f ; Mean nuclear charge\n", target->z2);
    printf("a2   = %10.3f ; Mean nuclear mass\n", target->a2);
    printf("iadj = %10.1f ; Logarithmic mean ionization potential [eV]\n", target->iadj);
    printf("rho  = %10.4e ; Density [g cm^-3]\n", target->rho);
    printf("pla  = %10.4e ; Plasma frequency [eV]\n", target->pla);
    printf("etad = %10.1f ; Ratio of density to density at STP for gasses (zero for everything else)\n", target->etad);
    printf("bind = %10.4e ; Total electronic binding energy [eV]\n", target->bind);
    printf("X0   = %10.4f ; Density effect turn-on value\n", target->X0);
    printf("X1   = %10.4f ; Density effect asymptotic bound\n", target->X1);
    printf("a    = %10.5f ; Density effect interpolation parameter\n", target->a);
    printf("m    = %10.4f ; Density effect interpolation parameter\n", target->m);
    printf("d0   = %10.2f ; Low energy density effect value (zero for everything but conductors)\n", target->d0);
    printf("\n");
    return;
}
