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
    return CRange::Tdata();
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
std::vector<CRange::Tdata> CRange::default_target(std::vector<CRange::Tdata> &extratargets)
{
    const std::string tnames[] = {"H", "He", "C", "N", "O", "Na", "Al", "Si", "P",
        "Ar", "Fe", "Ni", "Cu", "Ge", "Ag", "Ba", "Os", "Pt", "Au", "Pb", "U",
        "Air", "ArCO2", "BC-408", "BP-1", "CH2", "CO2", "CR-39", "CsI", "Halo",
        "Hosta", "ISM", "Kapton", "Kevlar", "Lexan", "LH2", "Mesh", "Mylar",
        "SiO2", "Teflon", "Water", "Unknown"};
    const double targets[][12] = {
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
///
/// \brief Allow printing.
///
std::ostream &operator<< (std::ostream &o, const CRange::Tdata &t) {
    t.print(&o);
    return o;
}
