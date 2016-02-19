///
/// \file crange.cpp
/// \brief Source code for the crange library functions.
///
/// This file contains all source code for the crange library.
///
//
// Include headers for variable and function declarations, &c..
//
#include "crange.h"
///
/// \brief Empty constructor.
///
CRange::Tdata::Tdata(void)
{
    _name = "Unknown";
    for (int i=0; i < Ndata; i++) data[i] = 0.0;
    init();
}
///
/// \brief Copy constructor.
///
/// Create a copy of a CRange::Tdata.
///
/// \param t The original CRange::Tdata object.
///
CRange::Tdata::Tdata( const CRange::Tdata &t )
{
    _name = t._name;
    for (int i=0; i < Ndata; i++) data[i] = t.data[i];
    _hash = t._hash;
}
///
/// \brief Standard constructor.
///
/// Initialize with a name and an array of data.
///
/// \param n Target name.
/// \param d Data array.
///
CRange::Tdata::Tdata( const std::string &n, const double d[] )
{
    _name = n;
    for (int i=0; i < Ndata; i++) data[i] = d[i];
    init();
}
///
/// \brief Standard constructor.
///
/// Initialize with a name and an array of data.
///
/// \param n Target name.
/// \param d Data array.
///
CRange::Tdata::Tdata( const char *n, const double d[] )
{
    _name = std::string(n);
    for (int i=0; i < Ndata; i++) data[i] = d[i];
    init();
}
///
/// \brief INI-based constructor.
///
/// Initialize a CRange::Tdata object with a section from an INI file.
///
/// \param n Target name; also a section of the INI file.
/// \param ini Pointer to an object representing the INI file.
///
CRange::Tdata::Tdata( const char *n, dictionary *ini )
{
    std::string namekey = std::string(n) + ":name";
    _name = iniparser_getstring(ini, namekey.c_str(), "Unknown");
    for (int i=0; i < Ndata; i++)
        data[i] = iniparser_getdouble(ini, (_name + ":" + dnames[i]).c_str(), 0.0);
    init();
}
///
/// \brief Consolidates TData initialization code.
///
/// This function contains the code that is common to all constructors.
///
void CRange::Tdata::init(void)
{
    _hash = 0.0;
    for (unsigned i=0; i < _name.length(); i++) _hash += (double)_name[i];
    for (int i=0; i < Ndata; i++) _hash += data[i];
}
///
/// \brief Assignment operator.
///
void CRange::Tdata::operator=( const CRange::Tdata &t )
{
    if (this == &t ) return;
    _name = t._name;
    for (int i=0; i < Ndata; i++) data[i] = t.data[i];
    _hash = t._hash;
}
///
/// \brief Print.
///
/// This serializes a CRange::Tdata object onto an output stream.
///
/// \param o Pointer to output stream.
///
void CRange::Tdata::print( std::ostream *o ) const
{
    *o << "[" << _name << "]" << std::endl;
    *o << "name = " << std::setw(10) << _name << " ; Target name" << std::endl;
    for (int i=0; i < Ndata; i++) {
        std::string n4 = dnames[i];
        for (int j = 4 - dnames[i].length(); j > 0; j--) n4 += " ";
        if (fixed[i]) {
            *o << n4 << " = " << std::setw(10) << std::fixed << std::setprecision(precision[i])
                << data[i] << " ; " << comment[i] << std::endl;
        } else {
            *o << n4 << " = " << std::setw(10) << std::scientific << std::setprecision(precision[i])
                << data[i] << " ; " << comment[i] << std::endl;
        }
    }
}
///
/// \brief Empty constructor.
///
CRange::RangeTable::RangeTable(void)
{
    z1 = a1 = 0.0;
    sswitch = 0;
    CRange::Tdata t;
    target = t;
    for (int i=0; i<MAXE; i++) range[i] = 0.0;
}
///
/// \brief Copy constructor.
///
CRange::RangeTable::RangeTable(const CRange::RangeTable &t)
{
    z1 = t.z1;
    a1 = t.a1;
    sswitch = t.sswitch;
    target = t.target;
    for (int i=0; i<MAXE; i++) range[i] = t.range[i];
}
///
/// \brief Assignment operator.
///
void CRange::RangeTable::operator=(const CRange::RangeTable &t)
{
    z1 = t.z1;
    a1 = t.a1;
    sswitch = t.sswitch;
    target = t.target;
    for (int i=0; i<MAXE; i++) range[i] = t.range[i];
}
///
/// \brief Standard constructor.
///
/// \param z Projectile charge.
/// \param a Projective mass number.
/// \param s Switches used.
/// \param t CRange::Tdata object representing the target material.
///
CRange::RangeTable::RangeTable(double z, double a, short s, CRange::Tdata &t)
{
    z1 = z;
    a1 = a;
    sswitch = s;
    target = t;
    int i = 0;
    while (CRange::energy_table(i) < 8.0) {
        double e1 = CRange::energy_table(i);
        range[i] = CRange::benton(e1, z1, a1, target);
        i++;
    }
    while (i < MAXE) {
        range[i] = range[i-1] + CRange::integrate_dedx(i, z1, a1, sswitch, target);
        i++;
    }
}
///
/// \brief Perform interpolation on an existing CRange::RangeTable.
///
/// \param e Initial projectile kinetic energy in A MeV.
///
/// \return Projectile range in g cm<sup>-2</sup>.
///
double CRange::RangeTable::interpolate_range( double e )
{
    if( e > CRange::energy_table(0)) {
        int i = 1;
        while ( e > CRange::energy_table(i) ) i++;
        return ( range[i-1] + ( e - CRange::energy_table(i-1) ) *
            (range[i]-range[i-1])/(CRange::energy_table(i)-CRange::energy_table(i-1)) );
    } else {
        return( e*range[0]/CRange::energy_table(0) );
    }
}
///
/// \brief Perform interpolation on an existing CRange::RangeTable.
///
/// \param e Initial projectile kinetic energy in A MeV.
/// \param r0 Thickness of material in g cm<sup>-2</sup>.
///
/// \return Final projectile kinetic energy after passing trough the material.
///
double CRange::RangeTable::interpolate_energy( double e, double r0 )
{
    double r;
    // Not sure why energy would ever be negative.
    if ( e > 0.0 ){
        //
        // Compute the expected total range at this energy.
        //
        double rr = interpolate_range(e);
        r = rr - r0;
        if (r < 0.0) return(0.0);
    } else {
        r = r0;
    }
    if (r > range[0]) {
        int i = 1;
        while( r > range[i] ) i++;
        return (energy_table(i-1) + (r - range[i-1]) * (energy_table(i)-energy_table(i-1)) / (range[i]-range[i-1]));
    } else {
        return CRange::energy_table(0) * r / range[0];
    }
}
///
/// \brief Computes effective projectile charge.
///
/// This is the modification of projectile charge due to electron
/// capture.  Hubert, Bimbot \& Gauvin, \cite art_fh2, give an
/// empirically determined function which depends on the target material.
/// This version is used if #SSWITCH_EC is set.
/// Two older versions, from Anthony \& Landford, \cite art_jma,
/// and Pierce \& Blann, \cite art_tep are also available.
///
/// \param z0 The bare projectile charge.
/// \param e1 The projectile kinetic energy in A MeV.
/// \param z2 The target mean nuclear charge.
/// \param sswitch The switch bit field.
///
/// \return The effective projectile charge.
///
/// \bug The Pierce \& Blann formula is not actually available; it is simply
/// commented out.
///
double CRange::effective_charge( double z0, double e1, double z2, short sswitch )
{
    double z1, capA, capB;
    double g = 1.0+e1/ATOMICMASSUNIT;
    double b2 = 1.0-1.0/(g*g);
    double b = sqrt(b2);
    double z23 = exp((2.0/3.0)*log(z0));
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
///
/// \brief Computes primary ionization.
///
/// This computes the primary ionization, the number of delta-rays produced
/// per unit length.  The formula is based on Bethe \cite art_hb,
/// as well as Fleischer <em>et al.</em>, \cite art_rlf3.
///
/// \param e1 The projectile kinetic energy in A MeV.
/// \param z0 The projectile charge.
/// \param I0 The binding energy of outermost electron in eV.
/// \param f0 The fraction of electrons in the outermost state.
/// \param K  A constant that depends on the target.
/// \param sswitch The switch bit field.
/// \param target A reference to a CRange::Tdata structure.
///
/// \return Number of delta-rays per unit length in units of g<sup>-1</sup> cm<sup>2</sup>.
///
/// \bug The parameters needed are not contained in the target table.
///
double CRange::djdx(double e1, double z0, double I0, double f0, double K, short sswitch, CRange::Tdata &target)
{
    double g = 1.0+e1/ATOMICMASSUNIT;
    double delt = ( sswitch & SSWITCH_ND ) ? delta(g,target) : olddelta(g,target);
    double b2 = 1.0-1.0/(g*g);
    // double b = sqrt(b2);
    double z2 = target.z2();
    double a2 = target.a2();
    double z1 = CRange::effective_charge(z0, e1, z2, sswitch);
    double f1 = 0.3070722*z1*z1*z2/(b2*a2);
    double f2 = log(2.0*ELECTRONMASS*b2*g*g/I0);
    double J = 0.5*f1*(f2-b2-delt+K)*(f0/I0);
    return J;
}
///
/// \brief Computes dE/dx.
///
/// This is the core of the whole package, the dE/dx calculator.  I have
/// based this largely on the work of Salamon, \cite tech_mhs.
/// Values of certain physical constants have been updated,
/// as well as some of the corrections to the basic stopping power formula.
///
/// If the restricted energy loss parameter \em rel0 is non-zero, dedx() computes
/// restricted energy loss instead.
///
/// The dE/dx calculator includes a number of effects that are controlled by
/// switches encoded in a bit field.  Below we describe each bit field and the
/// effect it controls.
///
///  - #SSWITCH_ND : Density effect version. If this bit is set (which it is
///    by default), a newer version of the density effect is used.  See
///    delta() and olddelta() for details.
///  - #SSWITCH_SH : Inner shell correction. The inner shell correction is
///    somewhat problematic.  It arises when the projectile velocity is
///    comparable to the velocity of inner shell electrons in the target medium.
///    This is discussed by Fano, \cite art_uf. The shell correction can be
///    included explicitly using this formula from Barkas \& Berger, \cite coll_whb.
///    Alternatively, the shell correction can be "hidden" in the logarithmic
///    mean ionization potential.  Much more work is required before this topic
///    can be fully understood.
///  - #SSWITCH_LE : Relativistic shell correction.  The Leung, or
///    relativistic shell correction is a small effect which is due to
///    relativistic inner shell electrons in very heavy targets.
///    See Leung, \cite art_ptl1, and Leung, \cite art_ptl2.  #SSWITCH_LE
///    has no effect unless #SSWITCH_SH is also turned on.
///  - The Lindhard-Sørensen effect (see lindhard()) is turned on by default.
///    The Bloch, Mott \& Ahlen effects are included for historical interestest.
///    Right now these can be turned on by uncommenting a particular section
///    of the code.
///  - #SSWITCH_KI : Ultrarelativistic kinematic correction.
///    This an estimate of the ultrarelativistic kinematic correction from
///    Ahlen, \cite art_spa2. It corrects to the
///    finite mass (as opposed to size) of the nucleus in relativistic
///    electron-nucleus collisions.
///  - #SSWITCH_RA : Radiative correction.
///    This is the radiative correction discussed in Ahlen, \cite art_spa2.
///    It arises from bremsstrahlung
///    of scattered electrons in ultrarelativistic collisions.  The
///    form here is that of Jankus, \cite art_vzj.
///    The parameter Q from that paper is here set equal to the geometric
///    mean between the the electron rest energy and \f$ 2 m_e c^2 \gamma \f$.
///  - #SSWITCH_PA : Slowing due to pair production. This value and the value for
///    the bremsstrahlung correction below are based on the work of
///    Sørensen, \cite coll_ahs.
///  - #SSWITCH_BR : Slowing due to projectile bremsstrahlung.  This version is
///    that of Sørensen, \cite coll_ahs, who has shown that this effect
///    is much smaller than the version suggested by Weaver \& Westphal,
///    \cite art_baw3..  This is due to their treatment of the projectile and
///    target nuclei as a point particles.  That version appeared in some
///    much older versions of this code, but has been replaced with
///    Sørensen's version.  We have not yet updated this code to reflect
///    Sørensen's more recent paper \cite art_ahs1.
///  - #SSWITCH_BA : Barkas effect.
///    This is the Barkas correction as calculated in Jackson \&
///    McCarthy, \cite art_jdj.  It is multiplied
///    by a factor of two to bring it into agreement with Lindhard, \cite art_jl1.
///    It is not, however, equal to the
///    results of Lindhard, and more work is needed to decide which, if any,
///    form is correct.  The recommended value seems to be the Jackson
///    \& McCarthy result multiplied by two.  Jackson \& McCarthy do not
///    have reliable values of \f$ F(V) \f$ for \f$ V < 0.8 \f$ .  For the purposes of the
///    computation, the cut-off is placed at \f$ V=1.0 \f$ .  I have followed the
///    convention of Salamon in having the Barkas correction multiply just
///    the "Bethe" portion of the stopping logarithm rather than the whole
///    stopping logarithm.  As there is considerable disagreement in the
///    literature about the application of correction, and as changing
///    the convention makes makes a difference of less than 1 A MeV even
///    in calculating the energy of stopping uranium, I have chosen to
///    leave it where it is.  Furthermore, I have found that a simple
///    power law \f$ V^{-2} \f$ is adequate to model Jackson \& McCarthy's function
///    for \f$ V > 1.0 \f$ , so I have used this instead of the numbers found by
///    reading off one of Jackson \& McCarthy's figures (these values are stored
///    in the array fva[10], but only the last value is used).
///
/// \param e1 The projectile kinetic energy in A MeV.
/// \param rel0 Restricted energy loss parameter in eV.
/// \param z0 The projectile charge.
/// \param a1 The projectile atomic number.
/// \param sswitch The switch bit field.
/// \param target A CRange::Tdata object.
///
/// \return dE/dx in units of A MeV g<sup>-1</sup> cm<sup>2</sup>
///
double CRange::dedx( double e1, double rel0, double z0, double a1, short sswitch, CRange::Tdata &target )
{
    const static double fva[10] = {0.33,0.078,0.03,0.014,0.0084,
                                   0.0053,0.0035,0.0025,0.0019,0.0014};
    double g = 1.0 + e1/ATOMICMASSUNIT;
    double delt = ( sswitch & SSWITCH_ND ) ? CRange::delta(g,target) : CRange::olddelta(g,target);
    double b2 = 1.0-1.0/(g*g);
    double b = sqrt(b2);
    double z2 = target.z2();
    double a2 = target.a2();
    double iadj = target.iadj();
    double z1 = CRange::effective_charge(z0, e1, z2, sswitch);
    double f1 = 0.3070722*z1*z1*z2/(b2*a1*a2);
    double f2 = log(2.0*ELECTRONMASS*b2/iadj);
    if( sswitch & SSWITCH_SH ){
        double etam2 = 1.0/(b*b*g*g);
        double cadj=1.0e-6*(iadj)*(iadj)*etam2*(0.422377
            +etam2*(0.0304043-etam2*0.00038106))
            +1.0e-9*(iadj)*(iadj)*(iadj)*etam2*(3.858019
            +etam2*(-0.1667989 + etam2*0.00157955));
        f2 -= cadj/(z2);
        if( sswitch & SSWITCH_LE ){
            f2 -= (5.0/3.0)*log(2.0*ELECTRONMASS*b2/iadj)*(1.0e+03*target.bind()/(z2*ELECTRONMASS))-(iadj*iadj/(4.0*ELECTRONMASS*ELECTRONMASS*b2));
        }
    }
    double f6 = 2.0*log(g)-b2;
    //
    // The Lindhard-Sorensen effect is now on by default.  The
    // Bloch-Mott-Ahlen effects are included for historical interest and
    // can be turned on by uncommenting the line after the next.
    //
    // double f3 = 0.0;
    double f3 = CRange::lindhard(z1,a1,b,sswitch); // Comment out this line if uncommenting the next.
    // double f3 = CRange::bma(z1,b);
    double f8 = ( sswitch & SSWITCH_KI ) ?
        0.5*(-log(1.0+2.0*((5.4858e-04)*g/a1)) - ((5.4858e-04)*g/a1)*b2/(g*g)) :
        0.0;
    double f9 = ( sswitch & SSWITCH_RA ) ?
        (ALPHA/M_PI)*b2*(6.0822 + log(2.0*g)*(log(2.0*g)*(2.4167 + 0.3333*log(2.0*g))-8.0314)) :
        0.0;
    double Spa = 0.0;
    if( sswitch & SSWITCH_PA ){
        double dpa = 1.0/sqrt(g);
        double ldpa=log(dpa);
        double l0 = log(2.0*g);
        // double Lpa0 = (19.0/9.0)*(log(g/4.0) - 11.0/6.0);
        double Lpa0s = (19.0/9.0)*log(183.0*exp(-1.0/3.0*log(z2))/(1.0 + 4.0*6.25470095193633*183.0*exp(-1.0/3.0*log(z2))/g));
        double Lpa1 = dpa*(4178.0/(81*M_PI*M_PI) - 21.0/27.0 - 248.0*l0/(27.0*M_PI*M_PI)
            +(28.0*l0/9.0 - 446.0/27.0)*2.0*ldpa/(M_PI*M_PI) + 14.0*4.0*ldpa*ldpa/(9.0*M_PI*M_PI));
        double Lpa = Lpa0s+Lpa1;
        Spa=4.08803936906434e-06*(z1*z1/a1)*(z2*z2/a2)*(1.0 + 1.0/z2)*g*Lpa;
    }
    double Sbr = 0.0;
    if ( sswitch & SSWITCH_BR ){
        double Bbr = log(1.0 + 2.0*g*0.179524783764566/(exp((1.0/3.0)*log(a1)) + exp((1.0/3.0)*log(a2)))/a1);
        Sbr = 5.21721169334564e-07*(z1*z1/a1)*(z1*z1/a1)*(z2*z2/a2)*g*Bbr;
    }
    double f4 = 1.0;
    if ( sswitch & SSWITCH_BA ) {
        double v = b*g/(ALPHA*sqrt(z2));
        if (v>1.0) {
            //
            // See the documentation for SSWITCH_BA for why this is commented
            // out.  Also note that the extrapolated power law is set to -2.0
            // instead of -2.5, as mentioned above.
            //
            // if (v < 9.0){
            //     i = 1;
            //     while (v >= (double)i) i++;
            //     fv = fva[i-1]+(v-(double)(i-1))*(fva[i]-fva[i-1]);
            // } else {
            //     fv = fva[9]*exp(-2.5*log(v/10.0));
            // }
            int i = 9;
            double fv = fva[i]*exp(-2.0*log(v/10.0));
            f4 = 1.0 + 2.0*z1*fv/(sqrt(z2));
        }
    }
    //
    // Compute restricted energy loss.  REL is activated by setting rel0 >0.
    //
    if (rel0 > 0.0) {
        double f7=log(2.0*ELECTRONMASS*b2*g*g/rel0)+b2*(rel0/(2.0*ELECTRONMASS*b2*g*g)-1.0);
        return f1*(f2*f4+f3+f6-delt/2.0 - 0.5*f7 +f8);
    } else {
        return f1*(f2*f4+f3+f6-(delt/2.0)+f8+f9) + Sbr + Spa;
    }
}
///
/// \brief Computes the density effect.
///
/// This function implements the density effect correction as formulated
/// in Sternheimer \& Peierls, \cite art_rms2 and as
/// extended in Sternheimer, Berger \& Seltzer, \cite art_rms1.
/// This version can distinguish between solids and gasses, and between
/// metals and insulators.  For conducting materials, there is a
/// low-energy density effect.
///
/// \param g Projectile Lorentz factor.
/// \param target A reference to a CRange::Tdata structure.
///
/// \return The value of the density effect.
///
double CRange::delta( double g, CRange::Tdata &target )
{
    double X0 = target.X0();
    double X1 = target.X1();
    double cbar = 2.0*log(target.iadj()/target.pla())+1.0;
    double b = sqrt(1.0 - 1.0/(g*g));
    double X = log10(b*g);
    double etad = target.etad();
    if (etad > 0) {
        cbar -= 2.303*log10(etad);
        X1 -= 0.5*log10(etad);
        X0 -= 0.5*log10(etad);
    }
    if (X < X0) {
        return (target.d0()*exp(4.6052*(X-X0)));
    } else if (X >= X0 && X < X1) {
        return (4.6052*X + exp(log(target.a()) + target.m()*log(X1-X)) - cbar);
    } else {
        return (4.6052*X - cbar);
    }
}
///
/// \brief Computes an obsolete version of the density effect.
///
/// This function implements the density effect correction as originally
/// formulated in Sternheimer \& Peierls, \cite art_rms2.
/// Although it is now obsolete, I have included it here for
/// compatibility with earlier codes.
///
/// \param g Projectile Lorentz factor.
/// \param target A reference to a CRange::Tdata structure.
///
/// \return The value of the density effect.
///
double CRange::olddelta( double g, CRange::Tdata &target )
{
    double y0,y1,dy3,a;
    if ( g < 1.8 ) return(0.0);
    double cbar = 2.0*log(target.iadj()/target.pla())+1.0;
    double b = sqrt(1.0 - 1.0/(g*g));
    double y = 2.0*log(b*g);
    double etad = target.etad();
    if( etad > 0 ){
        y+=log(etad);
        if(cbar>=12.25){
            y1=23.03;
            y0 = ( cbar >= 13.804 ) ? 1.502*cbar-11.52 : 9.212;
        } else {
            y1=18.42;
            if (cbar < 12.25) y0 = 9.212;
            if (cbar < 11.5) y0 = 8.751;
            if (cbar < 11.0) y0 = 8.291;
            if (cbar < 10.5) y0 = 7.830;
            if (cbar < 10.0) y0 = 7.370;
        }
    } else {
        if (target.iadj() >= 100.0){
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
///
/// \brief Computes the Bloch, Mott and Ahlen corrections.
///
/// This function computes the Mott correction of Ahlen, \cite art_spa1,
/// the Bloch correction of F. Bloch, \cite art_fb1,
/// and the Ahlen correction of Ahlen, \cite art_spa3.
/// All three of these corrections are
/// rendered obsolete by the Lindhard-Sørensen correction, and are
/// included here for historical interest and comparison with older
/// calculations.
///
/// \param z1 The projectile charge.
/// \param b The projectile velocity in units of the speed of light
/// (\em i.e. \f$ \beta = v/c \f$).
///
/// \return The sum of the Bloch, Mott and Ahlen corrections.
///
/// \note The variables lambda and theta0 are
/// free parameters in the Ahlen correction.  Theta0 also appears in
/// the Mott correction.  Here I have used Ahlen's recommended values,
/// lambda = 1, theta0 = 0.1.
/// An alternative formula, \f$ \theta_0 = \sqrt{\alpha/(\beta \gamma \lambda)} \f$ , is
/// suggested by Waddington, Freier \& Fixsen, \cite art_cjw1.
///
/// \warning The Mott correction has a severely
/// limited range of validity, especially for high charges.  It's so
/// bad it can render the calculation not just inaccurate, but
/// unphysical (dE/dx \< 0) below about 10 A MeV for uranium.  Ahlen
/// recommends turning the Mott correction off for \f$ Z/\beta > 100 \f$.
/// Here for \f$ Z/\beta > 100 \f$ the Mott correction is given the value at
/// \f$ Z/\beta = 100 \f$. This prescription is given by Waddington,
/// Freier \& Fixsen, \cite art_cjw1.
///
/// \bug Currently, this function is not called by anything.
///
double CRange::bma( double z1, double b )
{
    //
    // Compute a sum needed by the Bloch Correction
    //
    double y = z1*ALPHA/b;
    double y2 = y*y;
    int msum = (int)(10.0*y) + 1;
    double sumr = 0.0;
    for (int n=1; n<msum; n++) {
        double fn=(double)n;
        double fn2 = fn*fn;
        sumr += (1.0/(fn2+y2) - 1.0/fn2)/fn;
    }
    //
    // Compute the Bloch and Ahlen Corrections
    //
    double lambda = 1.0;
    double theta0 = 0.1;
    double f3 = -y2*(1.202+sumr) + CRange::relbloch(z1,b,lambda,theta0);
    //
    // The Mott term
    //
    Cdouble Cz1(0.5, -y);
    Cdouble Cz2(1.0, y);
    double cosx = cos(2.0 * (CRange::lngamma(Cz1).imag() +
                             CRange::lngamma(Cz2).imag()));
    double st = sin(theta0/2.0);
    double b2 = b*b;
    //
    // f5=0.5*z1a*(b*(1.725+(0.52-2.0*st*st)*M_PI*cosx)+
    //     z1a*(3.246-0.451*b2+z1a*(1.522*b+0.987/b+
    //     z1a*(4.569-0.494*b2-2.696/b2+
    //     z1a*(1.254*b+0.222/b
    //     -1.170/b2/b)))));
    //
    // if( (exp(9.0*log(y))/6.0) > fabs(f5/(f2*f4+f3+f5+f6-delt/2.0)) )
    if ( y > 100.0*ALPHA ) y = 100.0*ALPHA;
    double f5=0.5*b2*y*((1.725+(0.52-2.0*st)*M_PI*cosx)
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
///
/// \brief Compute the relativistic Bloch correction.
///
/// This is the relativistic Bloch (or Ahlen) correction of Ahlen,
/// \cite art_spa3.  The evaluation of this correction has
/// been enormously simplified by the use of fully complex arithmetic.
///
/// \param z12 The projectile charge.
/// \param b1 The projectile velocity in units of the speed of light
/// (\em i.e. \f$ \beta = v/c \f$).
/// \param lambda A free parameter, described in bma().
/// \param theta0 A free parameter, described in bma().
///
/// \return The value of the relativistic Bloch correction.
///
double CRange::relbloch( double z12, double b1, double lambda, double theta0 )
{
    const static double ge = 0.5772157;
    double nu = z12*ALPHA/b1;
    double g = 1.0/sqrt(1.0-b1*b1);
    double abgl=ALPHA/(b1*g*lambda);
    Cdouble Cnu1(1.0, nu);
    double sigma = CRange::lngamma(Cnu1).imag();
    Cdouble Ci(0.0,1.0);
    Cdouble Cnu(1.0, 2.0*nu);
    Cdouble Cth(theta0/2.0, 0.0);
    Cdouble Cabgl(abgl, 0.0);
    Cdouble Cmi(0.0,-1.0);
    Cdouble Cf1 = Cmi/Cnu;
    Cdouble Cf2 = std::pow(Cth, Cnu);
    Cdouble Cf3 = std::pow(Cabgl, Cnu);
    Cdouble Csigma2(0.0, 2.0*sigma);
    Cdouble Cf4 = std::exp(Csigma2);
    Cdouble Cone(1.0, 0.0);
    Cdouble Cf5 = Cone/Cnu;
    Cdouble Cf6 = std::log(Cth) - Cf5;
    Cdouble Cabglge(log(4.0/abgl)+ge-1.0, 0.0);
    Cdouble Cf7 = Cabglge + Cf5;
    Cdouble Cl1 = Cf1 * (2.0*Cf2 - Cf3*Cf4);
    Cdouble Cl2 = Cf1 * (2.0*Cf2*Cf6 + Cf3*Cf4*Cf7);
    double es = (M_PI*nu)*exp(M_PI*nu)/sinh(M_PI*nu);
    double bloch2 = (M_PI/2.0)*b1*b1*nu*es*(4.0*nu*log(2.0)*Cl1.real()
        + (nu*M_PI-1.0)*Cl1.imag() + 2.0*nu*Cl2.real());
    return bloch2;
}
///
/// \brief Compute the Lindhard-Sørensen correction.
///
/// This is the Lindhard-Sørensen correction including finite nuclear
/// size effects as described in Lindhard \& Sørensen, \cite art_jl2.
/// The defined variable #SSWITCH_NS will turn off the
/// nuclear size effect if it is set to zero.  For values of the Lorentz
/// factor above 10/R, where R is the nuclear size divided by the electron
/// Compton wavelength, the correction is set to its asymptotic value
/// which is described by Sørensen, \cite proc_ahs. This also avoids some
/// difficulties with the evaluation of the confluent hypergeometric function
/// (A. H. Sørensen, private communication).
///
/// \param zz The projectile charge.
/// \param aa The projectile atomic mass.
/// \param bb The projectile velocity in units of the speed of light
/// (\em i.e. \f$ \beta = v/c \f$).
/// \param sswitch The switch bit field.
///
/// \return The value of the Lindhard-Sørensen correction.
///
double CRange::lindhard( double zz, double aa, double bb, short sswitch )
{
    const static double compton=3.05573356675e-3; // 1.18 fm / Compton wavelength
    // Cdouble Cpi(M_PI, 0.0);
    // Cdouble Cone(1.0, 0.0);
    double a3 = exp(log(aa)/3.0);
    double eta = zz*ALPHA/bb;
    double gg = 1.0/sqrt(1.0 - bb*bb);
    double rho = a3*compton;
    double prh = bb*gg*rho;
    int n = 1;
    double sumterm = 0.0;
    double term1 = 0.0;
    double term3 = 1.0;
    double term2 = 0.0;
    // double pct = 1.0;
    //
    // Compute this section if the Lorentz factor is less than a certain
    // value OR if the Nuclear size switch is turned off (in which case
    // the Lorentz factor doesn't matter).
    //
    if ((gg < 10.0/rho) || !(sswitch & SSWITCH_NS)) {
        double dk[3];
        double dmk = 0.0;
        double dkm1 = 0.0;
        // while(fabs(pct) > 0.01) {
        while ( n < 100 ) {
            double k0 = (double)n;
            int max = n == 1 ? 3 : 2;
            for (int i=0; i<max; i++) {
                double k;
                if (i == 0) k = k0;
                if (i == 1) k = -k0-1.0;
                if (i == 2) k = -k0;
                double signk = k/fabs(k);
                double sk = sqrt(k*k - ALPHA*ALPHA*zz*zz);
                double l = (k>0) ? k : -k-1.0;
                Cdouble Cske(sk+1.0, eta);
                Cdouble Cketag(k, -eta/gg);
                Cdouble Cskmeta(sk, -eta);
                Cdouble Cexir = std::sqrt(Cketag/Cskmeta);
                Cdouble Cpiske(0.0,(M_PI/2.0)*(l-sk) - CRange::lngamma(Cske).imag());
                // Cdouble Cpiske(0.0,(M_PI/2.0)*(l-sk));
                Cdouble Cedr = Cexir*std::exp(Cpiske);
                double H = 0.0;
                Cdouble Ceds(0.0, 0.0);
                if ( sswitch & SSWITCH_NS ) {
                    Cdouble Cmske(-sk+1.0, eta);
                    Cdouble Cmskmeta(-sk, -eta);
                    Cdouble Cexis = std::sqrt(Cketag/Cmskmeta);
                    Cdouble Cpimske(0.0,(M_PI/2.0)*(l+sk) - CRange::lngamma(Cmske).imag());
                    Ceds = Cexis*std::exp(Cpimske);
                    Cdouble Caar = Cske;
                    Cdouble Caas = Cmske;
                    Cdouble Cbbr(2.0*sk + 1.0,0.0);
                    Cdouble Cbbs(-2.0*sk + 1.0,0.0);
                    Cdouble Czzr(0.0,2.0*prh);
                    Cdouble Cmprh(0.0,-prh);
                    Cdouble Clamr = Cexir*std::exp(Cmprh)*CRange::hyperg(Caar,Cbbr,Czzr);
                    Cdouble Clams = Cexis*std::exp(Cmprh)*CRange::hyperg(Caas,Cbbs,Czzr);
                    double grgs = Clamr.imag()/Clams.imag();
                    Cdouble Cgrgs = CRange::lngamma(Cbbs);
                    grgs *= exp( CRange::lngamma(Cske).real() -
                                 CRange::lngamma(Cmske).real() -
                                 CRange::lngamma(Cbbr).real() +
                                 Cgrgs.real() +
                                 2.0*sk*log(2.0*prh) );
                    if (cos(Cgrgs.imag())<1.0) grgs*= -1.0;
                    if (fabs(grgs) > 1.0e-9) {
                        double frgr = sqrt((gg-1.0)/(gg+1.0))*Clamr.real()/Clamr.imag();
                        double fsgs = sqrt((gg-1.0)/(gg+1.0))*Clams.real()/Clams.imag();
                        double gz = -1.0*signk*(rho*gg+1.5*ALPHA*zz);
                        double z1 = -1.0*signk*zz;
                        double b0 = 1.0;
                        double a0 = (1.0 + 2.0*fabs(k))*b0/(rho-gz);
                        double a1 = 0.5*(gz+rho)*b0;
                        double an = a1;
                        double anm1 = a0;
                        double bnm1 = b0;
                        double asum = a0;
                        double bsum = b0;
                        double nn = 1.0;
                        do {
                            double bn = ((rho-gz)*an + ALPHA*z1*anm1/2.0)/(2.0*nn+2.0*fabs(k)+1.0);
                            double anp1 = ((gz+rho)*bn - ALPHA*z1*bnm1/2.0)/(2.0*nn + 2.0);
                            asum += an;
                            bsum += bn;
                            nn += 1.0;
                            anm1 = an;
                            an = anp1;
                            bnm1 = bn;
                        } while(fabs(anm1/asum) > 1e-6 && fabs(bnm1/bsum) > 1e-6 );
                        double figi= (k>0) ? asum/bsum : bsum/asum;
                        H = ((frgr-figi)/(figi-fsgs))*grgs;
                    } else {
                        H = 0.0;
                    }
                }
                // double r = 1.0 + H*H + 2.0*H*(Cedr.real()*Ceds.real() + Cedr.imag()*Cedr.imag());
                dk[i] = std::arg(Cedr + Ceds*H);
            }
            if (n>1) dk[2] = dmk;
            double sdm2 = sin(dk[2]-dk[1]);
            term1 = k0*(k0+1.0)*sdm2*sdm2/(eta*eta*(2.0*k0 + 1.0));
            if (n>1) {
                double sd2 = sin(dk[0]-dkm1);
                term1 += k0*(k0-1.0)*sd2*sd2/(eta*eta*(2.0*k0 - 1.0));
            }
            double sdd = sin(dk[0]-dk[2]);
            term2 = k0*sdd*sdd/(eta*eta*(4.0*k0*k0 - 1.0));
            term3 = term1 - 1.0/k0;
            sumterm += term2 + term3;
            n += 1;
            dkm1 = dk[0];
            dmk = dk[1];
            // pct = (term2 + term3)/sumterm;
        }
    } else {
        sumterm = -log(prh) - 0.2; // Asymptotic value of the LS correction.
    }
    double lls = sumterm + 0.5*bb*bb;
    return lls;
}
///
/// \brief Compute a mathematical function related to bremsstrahlung.
///
/// This function is used in an obsolete version of projectile slowing
/// due to nuclear-nuclear bremsstrahlung.  It appears in Heitler's treatment
/// of bremsstrahlung, \cite book_wh, which was adapted by
/// Weaver \& Westphal, \cite art_baw3.
///
/// \param x The input parameter.
///
/// \return The value of the function.
///
/// \bug Currently, this function is unused.
///
double CRange::Fbrems( double x )
{
    int n = 1;
    double t = 1.0, s = 0.0;
    if (x == 1.0) {
        return M_PI*M_PI/12.0;
    } else if (x < 1.0) {
        while ((fabs(t) > 0.0001) || (n < 10)) {
            t *= -x;
            s += t/((double)(n*n));
            n++;
        }
        return -s;
    } else {
        while ((fabs(t) > 0.0001) || (n < 10)) {
            t *= -1.0/x;
            s += t/((double)(n*n));
            n++;
        }
        return (M_PI*M_PI/12.0 + 0.5*log(x)*log(x) + s);
    }
}
///
/// \brief Numerically integrate dedx() over one step.
///
/// To compute range, we integrate the reciprocal of dedx().  This
/// handles one step in the integration.
///
/// \param i Index into the range table.
/// \param z Projectile charge.
/// \param a Projective mass number.
/// \param s Switches used.
/// \param t CRange::Tdata object representing the target material.
///
/// \return Projectile range in g cm<sup>-2</sup>.
///
double CRange::integrate_dedx( int i, double z, double a, short s, CRange::Tdata &t)
{
    double rel = 0.0;
    double e0 = CRange::energy_table(i-1);
    double de2 = (CRange::energy_table(i) - e0)/2.0;
    double e1 = e0 + 1.33998104*de2;
    double dedx1 = CRange::dedx(e1,rel,z,a,s,t);
    double e2 = e0 + 1.86113631*de2;
    double dedx2 = CRange::dedx(e2,rel,z,a,s,t);
    double e3 = e0 + 0.13886369*de2;
    double dedx3 = CRange::dedx(e3,rel,z,a,s,t);
    double e4 = e0 + 0.66001869*de2;
    double dedx4 = CRange::dedx(e4,rel,z,a,s,t);
    double dr = de2 * (0.65214515/dedx1 +
                       0.34785485/dedx2 +
                       0.34785485/dedx3 +
                       0.65214515/dedx4);
    return dr;
}
///
/// \brief Computes total range given initial energy.
///
/// This function computes total range given initial energy.  The technique
/// is quite clever, in that if from one call to the next, the projectile
/// and target material parameters do not change, the calculation of
/// range is performed by table interpolation rather than direct integration.
/// The savings in calculation time can be enormous.  However, the
/// range of valid energies is limited by the size of the table.  The
/// function dE/dx is evaluated at most of the energies defined by the function
/// energy_table().  Results are stored in a vector of CRange::RangeTable objects.
///
/// \param e Initial projectile kinetic energy in A MeV.
/// \param z1 Projectile charge.
/// \param a1 Projectile atomic mass.
/// \param sswitch The switch bit field.
/// \param target A CRange::Tdata object.
/// \param rt A vector containing previously computed range tables.
///
/// \return Projectile range in g cm<sup>-2</sup>.
///
double CRange::range( double e, double z1, double a1, short sswitch, CRange::Tdata &target, std::vector<CRange::RangeTable> &rt )
{
    //
    // Search the range table for existing data
    //
    for (std::vector<CRange::RangeTable>::iterator it=rt.begin(); it != rt.end(); ++it) {
        if ((z1 == it->z1) && (a1 == it->a1) && (sswitch == it->sswitch) &&
            (target == it->target)) {
                return it->interpolate_range(e);
        }
    }
    //
    // If we didn't exit, we need to create a new table.
    //
    CRange::RangeTable table(z1, a1, sswitch, target);
    rt.push_back(table);
    return table.interpolate_range(e);
}
///
/// \brief Computes total range by direct integration of dE/dx.
///
/// This function computes total range by direct integration of the
/// dedx() function.  It does not create a range table or do table
/// interpolation.
///
/// \param e Initial energy in A MeV.
/// \param z1 Projectile charge.
/// \param a1 Projectile mass.
/// \param sswitch The switch bit field.
/// \param target A CRange::Tdata object.
///
/// \return Total range in g cm<sup>-2</sup>.
///
/// \bug The interpolation in the low-energy regime needs more attention.
/// \bug Currently, this function isn't called by anything.
///
double CRange::qrange( double e, double z1, double a1, short sswitch, CRange::Tdata &target )
{
    double ei = 8.0, en = 1.0;
    double ra[MAXE];
    if(e > ei){
        ra[0] = CRange::benton(ei,z1,a1,target);
        int i = 1;
        do {
            ra[i] = ra[i-1] + CRange::integrate_dedx(i, z1, a1, sswitch, target);
            i++;
        } while(en < e);
        return ra[i-1];
    } else if (e > 1.0 && e <= ei) {
        return CRange::benton(e,z1,a1,target);
    } else {
        return e*CRange::benton(en,z1,a1,target)/en;
    }
}
///
/// \brief Computes ranges at low energies.
///
/// This function is the result of empirical fits to very low energy
/// 1 A MeV \< E \< 8 A MeV ion ranges.  It follows the methods of
/// Barkas \& Berger, \cite coll_whb. A simplified
/// discussion, with a more complicated formula is given in
/// Benton \& Henke, \cite art_evb1.
/// As yet I know of no nicer way to deal with these low energies.
///
/// \param e Projectile kinetic energy in A MeV.
/// \param z1 Projectile charge.
/// \param a1 Projectile atomic mass.
/// \param target A CRange::Tdata object.
///
/// \return Projectile range in g cm<sup>-2</sup>.
///
/// \note The array join[4] demarcates three energy regions represented by
/// the three sets of coefficients in amn[3][4][4]. The demarcation is variable
/// in order to minimize discontinuities at the boundary.  The coefficients
/// in cjoin[2][7], which is used to initialize join[4], are inherited from
/// legacy code; I have not found them in the non-obscure literature.
/// Approximately, the three regions are \em E \< 1 A MeV, 1 \< \em E \< 7 A MeV
/// and \em E \> 7 A MeV.  I can find no reason why join[4] has four elements and
/// not two.
///
double CRange::benton( double e, double z1, double a1, CRange::Tdata &target )
{
    const static double amn[3][4][4] = {
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
    const static double czmn[4][4] = {
        {-6.000e-05,  5.252e-02,  1.285e-01,  0.000e+00},
        {-1.850e-03,  7.355e-02,  7.171e-02, -2.723e-02},
        {-7.930e-02,  3.323e-01, -1.234e-01,  1.530e-02},
        { 2.200e-01,  0.000e+00,  0.000e+00,  0.000e+00}
    };
    const static double cjoin[2][7] = {
        { 0.94,  20.19, -84.08,  132.98, -30.77, -102.29, 64.03},
        {12.62, -51.96, 199.14, -367.09, 327.06, -108.57,  0.00}
    };
    const double cr = PROTONMASS/ATOMICMASSUNIT;
    //
    // Compute join[4].
    //
    double join[4];
    for (int l = 0; l < 2; l++) {
        int m = 6;
        join[l] = cjoin[l][m];
        while (m > 0) join[l] = 0.001*target.iadj()*join[l] + cjoin[l][--m];
    }
    //
    // Compute prnglo[3].
    //
    double logt = log( e * cr );
    double logi = log( target.iadj() );
    double prnglo[3];
    for (int l = 0; l < 3; l++) {
        double loglambda=0.0;
        for (int m = 3; m >= 0; m-- ) {
            int n = 3;
            double term = amn[l][m][n];
            while (n > 0) term = term*logt + amn[l][m][--n];
            loglambda += term;
            loglambda *= logi;
        }
        loglambda /= logi;
        loglambda += log( (target.a2())/(target.z2()) );
        prnglo[l] = exp(loglambda);
        if (l == 1) prnglo[l] *= 1.0e-03;
    }
    //
    // Compute cz[4].
    //
    double g = 1.0 + e/ATOMICMASSUNIT;
    double b = sqrt(1.0 -1.0/(g*g));
    double x = 137.0*b/z1;
    double cz[4];
    for (int m = 0; m < 4; m++) {
        cz[m]=0.0;
        for (int n = 3; n >= 0; n--) {
            cz[m] += czmn[m][n];
            cz[m] *= x;
        }
        cz[m] /= x;
    }
    int l = 1;
    if ( e < join[0]) l = 0;
    if ( e > join[1]) l = 2;
    int n;
    if (x <= 0.2) {
        n = 0;
    } else if ( x > 0.2 && x <= 2.0 ) {
        n = 1;
    } else if ( x > 2.0 && x <= 3.0 ) {
        n = 2;
    } else { // if ( x > 3.0 ) {
        n = 3;
    }
    double bzz=(31.8+3.86*exp((5.0/8.0)*logi))
        *( (target.a2())/(target.z2()) )*1.0e-06*exp((8.0/3.0)*log( z1 ));
    return(( (a1/cr)/(z1*z1) )*(prnglo[l] + bzz*cz[n]));
}
///
/// \brief Extract energies from range tables.
///
/// This function extracts energies from a range table by table interpolation.
///
/// \param e Projectile kinetic energy [A MeV].
/// \param r0 Range [g cm<sup>-2</sup>].
/// \param z1 Projectile charge.
/// \param a1 Projectile atomic mass.
/// \param sswitch The switch bit field.
/// \param target A CRange::Tdata object.
/// \param rt A vector containing previously computed range tables.
///
/// \return The final energy of the projectile.
///
double CRange::renergy( double e, double r0, double z1, double a1, short sswitch, CRange::Tdata &target, std::vector<CRange::RangeTable> &rt )
{
    //
    // Search the range table for existing data
    //
    for (std::vector<CRange::RangeTable>::iterator it=rt.begin(); it != rt.end(); ++it) {
        if ((z1 == it->z1) && (a1 == it->a1) && (sswitch == it->sswitch) &&
            (target == it->target)) {
                return it->interpolate_energy(e, r0);
        }
    }
    //
    // If we didn't exit, we need to create a new table.
    //
    CRange::RangeTable table(z1, a1, sswitch, target);
    rt.push_back(table);
    return table.interpolate_energy(e, r0);
}
///
/// \brief Print usage message.
///
/// Prints a usage message on STDERR.
///
/// \param executable Name of the program, usually argv[0].
///
void CRange::usage( char *executable )
{
    std::cerr << "usage: " << executable << " [-c COMMAND] [-h] [-l] [-o FILE] [-s switch.ini] [-t target.ini] [-V] <task file>" << std::endl;
    std::cerr << "       -c COMMAND    = Execute this one-line command instead of reading it from a file." << std::endl;
    std::cerr << "       -h            = Print this help message and exit." << std::endl;
    std::cerr << "       -l            = Print the built-in target table and exit." << std::endl;
    std::cerr << "       -o FILE       = Write to this file instead of standard output." << std::endl;
    std::cerr << "       -s switch.ini = Override the default switch values by reading this file." << std::endl;
    std::cerr << "       -t target.ini = Override the default target values by reading this file." << std::endl;
    std::cerr << "       -V            = Print the version string and exit." << std::endl;
    std::cerr << "       <task file>   = A file containing a list of tasks for crange.  Required unless a command is specified with -c." << std::endl;
}
///
/// \brief Print version message.
///
/// Prints a version message on STDERR.
///
/// \param executable Name of the program, usually argv[0].
///
void CRange::version(char *executable)
{
    std::cerr << executable << " version "
              << CRange_VERSION_MAJOR << "."
              << CRange_VERSION_MINOR << "."
              << CRange_VERSION_PATCH << "." << std::endl;
}
///
/// \brief Parses and executes the task list.
///
/// This utility function steps through the range, energy and dE/dx tasks
/// specified in the input data file.
/// The tasks are denoted by a single letter:
///     - \em r compute ranges
///     - \em e compute energies
///     - \em d compute dE/dx
///     - \em j compute dJ/dx (primary ionization)
///
/// The task letter should be followed by the energy (or range) at which to
/// compute range (or energy), the charge and mass of the particle, and the
/// name of the target material.  Names of target materials can be found in the
/// target.ini file.  Target material names should contain no whitespace.
///
/// \param commands A set of commands to execute
/// \param sswitch The switch bit field.
/// \param targets A vector of CRange::Tdata.
///
/// \return A vector of numerical results, ready for output.
///
/// \bug The primary ionization parameters are currently hard-coded.
/// \bug The numerical format and precision are currently hard-coded.
///
std::vector<std::string> CRange::process( std::vector<std::string> &commands, short sswitch, std::vector<CRange::Tdata> &targets )
{
    std::vector<std::string> results;
    if (commands.size() == 0) return results;
    std::vector<CRange::RangeTable> rt;
    for (std::vector<std::string>::iterator it=commands.begin(); it != commands.end(); ++it) {
        std::istringstream c(*it);
        std::ostringstream r;
        std::string task, targname;
        double red1, red2, z1, a1;
        double out = -9999.0;
        c >> task >> red1 >> red2 >> z1 >> a1 >> targname;
        //
        // c.good() will be false, even on a good read, because a good read
        // will also reach the end of input, so we test eof instead.
        //
        if (c.eof()) {
            CRange::Tdata target = CRange::find_target(targname, targets);
            if (target.name() == "Unknown") {
                std::cerr << "Invalid target detected in command: " << *it << std::endl;
            } else {
                if (task == "r") {
                    out = CRange::range(red1,z1,a1,sswitch,target,rt);
                } else if (task == "e") {
                    out = CRange::renergy(red1,red2,z1,a1,sswitch,target,rt);
                } else if (task == "d") {
                    out = CRange::dedx(red1,red2,z1,a1,sswitch,target);
                } else if (task == "j") {
                    out = CRange::djdx(red1, z1, 2.0, 0.05, 3.04, sswitch, target);
                } else {
                    std::cerr << "Invalid task detected in command: " << *it << std::endl;
                }
            }
        } else {
            std::cerr << "Invalid command detected: " << *it << std::endl;
        }
        r << std::setprecision(6) << std::fixed << out;
        // r << std::setprecision(5) << std::scientific << out;
        results.push_back(r.str());
    }
    return results;
}
///
/// \brief Initializes the value of of the switch bit field.
///
/// This utility reads an INI-type file and sets the switch bit field
/// accordingly.
///
/// \param switchfile The name of an INI-type file containing switch configuration.
///
/// \return The switch bit field.
///
/// \warning If the file is not found, this function will only
/// return the default value #SSWITCH_DEFAULT.
///
short CRange::init_switch( const std::string &switchfile )
{
    short sswitch = SSWITCH_DEFAULT;
    if (switchfile.length() == 0) return sswitch;
    if (access(switchfile.c_str(), R_OK) == 0) {
        sswitch = 0;
        dictionary *d = iniparser_load( switchfile.c_str() );
        if (iniparser_getboolean(d,"switch:barkas",0)) sswitch |= SSWITCH_BA;
        if (iniparser_getboolean(d,"switch:shell",0))  sswitch |= SSWITCH_SH;
        if (iniparser_getboolean(d,"switch:leung",0))  sswitch |= SSWITCH_LE;
        if (iniparser_getboolean(d,"switch:new delta",1))  sswitch |= SSWITCH_ND; // True by default!
        if (iniparser_getboolean(d,"switch:new electron capture",0))  sswitch |= SSWITCH_EC;
        if (iniparser_getboolean(d,"switch:finite nuclear size",1))  sswitch |= SSWITCH_NS; // True by default!
        if (iniparser_getboolean(d,"switch:kinematic",0))  sswitch |= SSWITCH_KI;
        if (iniparser_getboolean(d,"switch:radiative",0))  sswitch |= SSWITCH_RA;
        if (iniparser_getboolean(d,"switch:pair",0))  sswitch |= SSWITCH_PA;
        if (iniparser_getboolean(d,"switch:bremsstrahlung",0))  sswitch |= SSWITCH_BR;
        iniparser_freedict(d);
    } else {
        std::cerr << "Could not read switch file: " << switchfile
                  << ". Using default crange switches." << std::endl;
    }
    return sswitch;
}
///
/// \brief Initializes the value of of the switch bit field.
///
/// This utility reads an INI-type file and sets the switch bit field
/// accordingly.
///
/// \param switchfile The name of an INI-type file containing switch configuration.
///
/// \return The switch bit field.
///
/// \warning  If the file is not found, this function will only
/// return the default value #SSWITCH_DEFAULT.
///
short CRange::init_switch( const char *switchfile )
{
    const std::string switchstring = switchfile;
    return CRange::init_switch(switchstring);
}
///
/// \brief Read optional target data file.
///
/// This utility reads an INI-type file and returns a vector of pointers to
/// CRange::Tdata objects.
///
/// \param targetfile the name of an INI-type file containing target data.
///
/// \return A vector of CRange::Tdata objects.
///
/// \warning If the file is not found, the vector will be empty.
///
std::vector<CRange::Tdata> CRange::init_target( const std::string &targetfile )
{
    std::vector<CRange::Tdata> target_list;
    if (targetfile.length() > 0) {
        if (access(targetfile.c_str(), R_OK) == 0) {
            dictionary *d = iniparser_load(targetfile.c_str());
            int nsec = iniparser_getnsec(d);
            for (int i = 0; i < nsec; i++) {
                const char *sec = iniparser_getsecname(d, i);
                CRange::Tdata t(sec, d);
                target_list.push_back(t);
            }
            iniparser_freedict(d);
        } else {
            std::cerr << "Could not read target file: " << targetfile
                      << ". Using default target data." << std::endl;
        }
    }
    target_list = CRange::default_target(target_list);
    return target_list;
}
///
/// \brief Read optional target data file.
///
/// This utility reads an INI-type file and returns a vector of pointers to
/// CRange::Tdata objects.
///
/// \param targetfile the name of an INI-type file containing target data.
///
/// \return A vector of CRange::Tdata objects.
///
/// \warning If the iniparser library is not found, this function will only
/// return a NULL pointer.
///
std::vector<CRange::Tdata> CRange::init_target( const char *targetfile )
{
    const std::string targetstring = targetfile;
    return CRange::init_target(targetstring);
}
///
/// \brief Find a target by name.
///
CRange::Tdata CRange::find_target(const char *name, std::vector<CRange::Tdata> &targets)
{
    const std::string strname = name;
    return CRange::find_target(strname, targets);
}
///
/// \brief Find a target by name.
///
CRange::Tdata CRange::find_target(const std::string &name, std::vector<CRange::Tdata> &targets)
{
    for (std::vector<CRange::Tdata>::iterator it=targets.begin(); it != targets.end(); ++it) {
        if (it->name() == name) {
            return *it;
        }
    }
    CRange::Tdata t;
    return t;
}
///
/// \brief Load the default target data.
///
/// This function returns a vector of CRange::Tdata containing the default
/// target data. The default data may be added to or overridden by supplying an
/// INI-type file on the command line.  This is accomplished by loading the
/// extra target data into the vector first, then loading the default data.
///
/// \param extratargets The targets already loaded.
///
/// \return A vector of CRange::Tdata.
///
std::vector<CRange::Tdata>& CRange::default_target(std::vector<CRange::Tdata> &extratargets)
{
    const static std::string tnames[] = {"H", "He", "C", "N", "O", "Na", "Al", "Si", "P",
        "Ar", "Fe", "Ni", "Cu", "Ge", "Ag", "Ba", "Os", "Pt", "Au", "Pb", "U",
        "Air", "ArCO2", "BC-408", "BP-1", "CH2", "CO2", "CR-39", "CsI", "Halo",
        "Hosta", "ISM", "Kapton", "Kevlar", "Lexan", "LH2", "Mesh", "Mylar",
        "SiO2", "Teflon", "Water", "Unknown"};
    const static double targets[][12] = {
        //    z2,      a2,  iadj,        rho,        pla,etad,       bind,      X0,     X1,       a,      m,   d0
        {  1.000,   1.008,  19.2, 8.3748e-05, 2.6300e-01, 1.0, 1.3606e-02,  1.8639, 3.2718, 0.14092, 5.7273, 0.00 },
        {  2.000,   4.003,  41.8, 1.6632e-04, 2.6300e-01, 1.0, 7.7872e-02,  2.2017, 3.6122, 0.13443, 5.8347, 0.00 },
        {  6.000,  12.011,  78.0, 2.0000e+00, 2.8803e+01, 0.0, 1.0251e+00, -0.0351, 2.4860, 0.20240, 3.0036, 0.10 },
        {  7.000,  14.007,  82.0, 1.1653e-03, 6.9500e-01, 1.0, 1.4782e+00,  1.7378, 4.1323, 0.15349, 3.2125, 0.00 },
        {  8.000,  15.999,  95.0, 1.3315e-03, 7.4400e-01, 0.0, 2.0359e+00,  1.7541, 4.3213, 0.11778, 3.2913, 0.00 },
        { 11.000,  22.990, 149.0, 9.7100e-01, 1.9641e+01, 0.0, 4.4097e+00,  0.2880, 3.1962, 0.07772, 3.6452, 0.08 },
        { 13.000,  26.980, 166.0, 2.6989e+00, 3.2860e+01, 0.0, 6.5929e+00,  0.1708, 3.0127, 0.08024, 3.6345, 0.12 },
        { 14.000,  28.090, 173.0, 2.3300e+00, 3.1055e+01, 0.0, 7.8751e+00,  0.2014, 2.8715, 0.14921, 3.2546, 0.14 },
        { 15.000,  30.974, 173.0, 2.2000e+00, 2.9743e+01, 0.0, 9.2905e+00,  0.1696, 2.7815, 0.23610, 2.9158, 0.14 },
        { 18.000,  39.948, 188.0, 1.6620e-03, 7.8900e-01, 1.0, 1.4382e+01,  1.7635, 4.4855, 0.19714, 2.9618, 0.00 },
        { 26.000,  55.847, 286.0, 7.8740e+00, 5.5172e+01, 0.0, 3.4582e+01, -0.0012, 3.1531, 0.14680, 2.9632, 0.12 },
        { 28.000,  58.690, 311.0, 8.9020e+00, 5.9385e+01, 0.0, 4.1324e+01, -0.0566, 3.1851, 0.16496, 2.8430, 0.10 },
        { 29.000,  63.550, 323.0, 8.9600e+00, 5.8270e+01, 0.0, 4.4973e+01, -0.0254, 3.2792, 0.14339, 2.9044, 0.08 },
        { 32.000,  72.610, 350.0, 5.3230e+00, 4.4141e+01, 0.0, 5.7047e+01,  0.3376, 3.6096, 0.07188, 3.3306, 0.14 },
        { 47.000, 107.870, 470.0, 1.0500e+01, 6.1635e+01, 0.0, 1.4451e+02,  0.0657, 3.1074, 0.24585, 2.6899, 0.14 },
        { 56.000, 137.330, 491.0, 3.5000e+00, 3.4425e+01, 0.0, 2.2118e+02,  0.4190, 3.4547, 0.18268, 2.8906, 0.14 },
        { 76.000, 190.200, 746.0, 2.2570e+01, 8.6537e+01, 0.0, 4.6939e+02,  0.0891, 3.5414, 0.12751, 2.9608, 0.10 },
        { 78.000, 195.090, 790.0, 2.1450e+01, 8.4389e+01, 0.0, 5.0101e+02,  0.1484, 3.6212, 0.11128, 3.0417, 0.12 },
        { 79.000, 196.970, 790.0, 1.9320e+01, 8.0215e+01, 0.0, 5.1732e+02,  0.2021, 3.6979, 0.09756, 3.1101, 0.14 },
        { 82.000, 207.200, 823.0, 1.1350e+01, 6.1072e+01, 0.0, 5.6834e+02,  0.3776, 3.8073, 0.09359, 3.1608, 0.14 },
        { 92.000, 238.030, 890.0, 1.5370e+01, 7.7986e+01, 0.0, 7.6220e+02,  0.2260, 3.3721, 0.19677, 2.8171, 0.14 },
        {  7.312,  14.667,  85.4, 1.2048e-03, 7.0700e-01, 1.0, 1.7150e+00,  1.7418, 4.2759, 0.10914, 3.3994, 0.00 },
        { 15.867,  34.892, 174.7, 1.7000e-03, 8.0120e-01, 1.0, 1.7150e+00,  1.7418, 4.2759, 0.10914, 3.3994, 0.00 },
        {  3.381,   6.248,  62.8, 1.0320e+00, 2.1534e+01, 1.0, 4.9530e-01,  0.1769, 2.6747, 0.11442, 3.3762, 0.00 },
        { 14.795,  32.575, 242.5, 3.0000e+00, 3.3636e+01, 0.0, 2.7489e+01,  0.0843, 3.6297, 0.06445, 3.3655, 0.00 },
        {  2.667,   4.676,  57.4, 9.4000e-01, 2.1099e+01, 0.0, 3.5078e-01,  0.1370, 2.5177, 0.12108, 3.4292, 0.00 },
        {  7.333,  14.670,  85.0, 1.8421e-03, 8.7400e-01, 1.0, 1.6990e+00,  1.6294, 4.1825, 0.11768, 3.3227, 0.00 },
        {  3.946,   7.413,  73.1, 1.4000e+00, 2.4780e+01, 0.0, 7.2426e-01,  0.1562, 2.6507, 0.12679, 3.3076, 0.00 },
        { 54.000, 129.905, 553.1, 4.5100e+00, 3.9455e+01, 0.0, 2.0258e+02,  0.0395, 3.3353, 0.25381, 2.6657, 0.00 },
        {  1.077,   1.238,  19.2, 4.0880e-24, 5.4342e-11, 1.0, 1.8554e-02,  2.0000, 3.5000, 0.13500, 5.7500, 0.00 },
        {  4.250,   8.091,  72.3, 1.4000e+00, 2.4722e-01, 0.0, 7.7212e-01,  0.1606, 2.6255, 0.12860, 3.3288, 0.00 },
        {  1.077,   1.238,  19.2, 2.0440e-24, 3.8426e-11, 1.0, 1.8554e-02,  2.0000, 3.5000, 0.13500, 5.7500, 0.00 },
        {  5.026,   9.803,  75.9, 1.4200e+00, 2.4586e+01, 0.0, 9.1859e-01,  0.1509, 2.5631, 0.15972, 3.1912, 0.00 },
        {  4.000,   7.567,  71.7, 1.4500e+00, 2.5229e+01, 0.0, 9.1859e-01,  0.1509, 2.5631, 0.15972, 3.1912, 0.00 },
        {  4.061,   7.706,  73.1, 1.2040e+00, 2.2915e+01, 0.0, 6.8789e-01,  0.1606, 2.6255, 0.12860, 3.3288, 0.00 },
        {  1.000,   1.008,  21.8, 6.0000e-02, 7.0310e+00, 0.0, 1.3606e-02,  0.4759, 1.9215, 0.13483, 5.6249, 0.00 },
        { 21.400,  44.200, 223.0, 2.1500e+00, 2.9400e+01, 0.0, 3.4582e+01, -0.0012, 3.1531, 0.14680, 2.9632, 0.12 },
        {  4.456,   8.735,  78.7, 1.4000e+00, 2.4595e+01, 0.0, 8.4108e-01,  0.1562, 2.6507, 0.12679, 3.3067, 0.00 },
        { 10.000,  20.029, 139.2, 2.3200e+00, 3.1014e+01, 0.0, 3.9822e+00,  0.1385, 3.0025, 0.08408, 3.5064, 0.00 },
        {  8.000,  16.669,  99.1, 2.2000e+00, 2.9609e+01, 0.0, 2.1465e+00,  0.1648, 2.7404, 0.10606, 3.4046, 0.00 },
        {  3.333,   6.005,  75.0, 1.0000e+00, 2.1469e+01, 0.0, 6.8770e-01,  0.2400, 2.8004, 0.09116, 3.4773, 0.00 },
        // THIS MUST BE THE LAST STRUCTURE DEFINITION!
        {  0.000,   0.000,   0.0, 0.0000e+00, 0.0000e+00, 0.0, 0.0000e+00,  0.0000, 0.0000, 0.00000, 0.0000, 0.00 }
    };
    int k=0;
    do {
        CRange::Tdata t(tnames[k], targets[k]);
        extratargets.push_back(t);
        k++;
    } while (tnames[k] != "Unknown");
    CRange::Tdata t(tnames[k], targets[k]); // Add "Unknown" to the end.
    extratargets.push_back(t);
    return extratargets;
}
