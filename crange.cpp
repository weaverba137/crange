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
}
///
/// \brief Copy constructor.
///
CRange::Tdata::Tdata( const CRange::Tdata &t )
{
    _name = t._name;
    for (int i=0; i < Ndata; i++) data[i] = t.data[i];
}
///
/// \brief Standard constructor.
///
CRange::Tdata::Tdata( const std::string &n, const double d[] )
{
    _name = n;
    for (int i=0; i < Ndata; i++) data[i] = d[i];
}
///
/// \brief Slightly different constructor.
///
CRange::Tdata::Tdata( const char *n, const double d[] )
{
    _name = std::string(n);
    for (int i=0; i < Ndata; i++) data[i] = d[i];
}
///
/// \brief INI based constructor.
///
CRange::Tdata::Tdata( const char *n, dictionary *ini )
{
    std::string namekey = std::string(n) + ":name";
    _name = iniparser_getstring(ini, namekey.c_str(), "Unknown");
    for (int i=0; i < Ndata; i++) {
        data[i] = iniparser_getdouble(ini, (_name + ":" + dnames[i]).c_str(), 0.0);
    }
}
///
/// \brief Print.
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
    double b = sqrt(b2);
    double z2 = target.z2();
    double a2 = target.a2();
    double z1 = CRange::effective_charge(z0, e1, z2, sswitch);
    double f1 = 0.3070722*z1*z1*z2/(b2*a2);
    double f2 = log(2.0*ELECTRONMASS*b2*g*g/I0);
    double J = 0.5*f1*(f2-b2-delt+K)*(f0/I0);
    return J;
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
              << crange_VERSION_MAJOR << "."
              << crange_VERSION_MINOR << "."
              << crange_VERSION_PATCH << "." << std::endl;
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
/// \warning If the iniparser library is not found, this function will only
/// return the default value #SSWITCH_DEFAULT.
///
short CRange::init_switch( const char *switchfile )
{
    short sswitch = SSWITCH_DEFAULT;
    if (access(switchfile,R_OK) == 0) {
        sswitch = 0;
        dictionary *d = iniparser_load( switchfile );
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
/// \warning If the iniparser library is not found, this function will only
/// return the default value #SSWITCH_DEFAULT.
///
short CRange::init_switch( const std::string &switchfile )
{
    return CRange::init_switch(switchfile.c_str());
}
///
/// \brief Read optional target data file.
///
/// This utility reads an INI-type file and returns an array of pointers to
/// CRange::Tdata structures.
///
/// \param targetfile the name of an INI-type file containing target data.
///
/// \return A pointer to an array of CRange::Tdata structures.  This pointer must
/// be free()d!
///
/// \warning If the iniparser library is not found, this function will only
/// return a NULL pointer.
///
std::vector<CRange::Tdata> CRange::init_target( const char *targetfile )
{
    std::vector<CRange::Tdata> target_list;
    dictionary *ini = iniparser_load(targetfile);
    int nsec = iniparser_getnsec(ini);
    for (int i = 0; i < nsec; i++) {
        const char *sec = iniparser_getsecname(ini, i);
        CRange::Tdata t(sec, ini);
        target_list.push_back(t);
    }
    return target_list;
}
///
/// \brief Read optional target data file.
///
/// This utility reads an INI-type file and returns an array of pointers to
/// CRange::Tdata structures.
///
/// \param targetfile the name of an INI-type file containing target data.
///
/// \return A pointer to an array of CRange::Tdata structures.  This pointer must
/// be free()d!
///
/// \warning If the iniparser library is not found, this function will only
/// return a NULL pointer.
///
std::vector<CRange::Tdata> CRange::init_target( const std::string &targetfile )
{
    return CRange::init_target(targetfile.c_str());
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
    return CRange::Tdata();
}
///
/// \brief Allow printing.
///
std::ostream &operator<< (std::ostream &o, const CRange::Tdata &t) {
    t.print(&o);
    return o;
}
