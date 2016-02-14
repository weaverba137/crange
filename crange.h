///
/// \file crange.h
/// \brief Header file for crange.
///
/// This header file collects all the other header files needed to compile
/// crange, as well as all defines, function declarations, etc.
///
//
// Standard includes here to keep crange.c clean
//
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <unistd.h> // to get getopt()
#include <cerrno>
#include <ctime>
//
// Include config.h
//
#include <config.h>
//
// Include header for parsing ini files.
//
#ifdef HAVE_INIPARSER_H
#include <iniparser.h>
#endif
//
// These values define the range-energy tables
//
#define LOGTENEMIN 0.0 ///< \f$ \log_{10} E_{\mathrm{min}} \f$ Minimum energy in units of A MeV.
#define LOGTENEMAX 6.0 ///< \f$ \log_{10} E_{\mathrm{max}} \f$ Maximum energy in units of A MeV.
#define MAXE 200 ///< The number of energies in the range-energy tables.
#define MAXAB 50 ///< The number of range tables.  Arbitrary, but should be larger than the number of targets.
//
// M_PI, M_PI_2 and M_LN10 should be defined in math.h.
//
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 ///< The value of pi, in case it is not defined in math.h.
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169163975144 ///< The value of pi/2, in case it is not defined in math.h.
#endif
#ifndef M_LN10
#define M_LN10 2.30258509299404568402 ///< The value of \f$ \ln 10 \f$ in case it is not defined in math.h.
#endif
//
// Physical constants.
//
#ifdef ALPHA
#undef ALPHA
#endif
#ifndef ALPHA
#define ALPHA 7.29735301383e-3 ///< The fine structure constant.
#endif
#ifdef ATOMICMASSUNIT
#undef ATOMICMASSUNIT
#endif
#ifndef ATOMICMASSUNIT
#define ATOMICMASSUNIT 931.4943 ///< 1 amu in units of MeV/c<sup>2</sup>.
#endif
#ifdef PROTONMASS
#undef PROTONMASS
#endif
#ifndef PROTONMASS
#define PROTONMASS 938.2723 ///< Proton mass in units of MeV/c<sup>2</sup>.
#endif
#ifdef ELECTRONMASS
#undef ELECTRONMASS
#endif
#ifndef ELECTRONMASS
#define ELECTRONMASS 0.511003e+6 ///< Electron mass in units of eV/c<sup>2</sup>.
#endif
//
// These define the values in the switch bit field.
//
#define SSWITCH_BA  0x001 ///< Barkas effect bit.
#define SSWITCH_SH  0x002 ///< Shell effect bit.
#define SSWITCH_LE  0x004 ///< Leung effect bit.
#define SSWITCH_ND  0x008 ///< Density effect switcher bit.
#define SSWITCH_EC  0x010 ///< Electron capture switcher bit.
#define SSWITCH_NS  0x020 ///< Finite Nuclear Size effect bit.
#define SSWITCH_KI  0x040 ///< Kinetic effect bit.
#define SSWITCH_RA  0x080 ///< Radiative correction effect bit.
#define SSWITCH_PA  0x100 ///< Pair Production energy loss bit.
#define SSWITCH_BR  0x200 ///< Projectile Bremsstrahlung bit.
#define SSWITCH_DEFAULT (SSWITCH_ND | SSWITCH_NS) ///< Default bits set density effect and finite nuclear size.
//
// Probably can switch to using std::string.
//
#define NAMEWIDTH 8 ///< The maximum number of characters in a target name.
//
// Time to define the namespace.
//
namespace CRange
{
///
/// \brief Structure containing target data.
///
/// This structure contains all the data related to target materials.
///
struct TDATA {
    ///
    /// \name Material name
    /// @{
    ///
    char name[NAMEWIDTH+1]; ///< The name of the material.
    /// @}
    ///
    /// \name General parameters
    /// @{
    ///
    double z2;   ///< The mean nuclear charge.
    double a2;   ///< The mean atomic mass number.
    double iadj; ///< The logarithmic mean ionization potential [eV].
    double rho;  ///< The density [g cm<sup>-3</sup>].
    double pla;  ///< The plasma frequency [eV].
    double etad; ///< The ratio of density to density at STP for gaseous targets. Should be set to zero for non-gaseous materials.
    double bind; ///< The total electronic binding energy [keV].
    //// @}
    ///
    /// \name Density effect parameters
    /// @{
    ///
    double X0; ///< Value of \f$ \log_{10} \beta\gamma \f$ at which the density effect turns on.
    double X1; ///< Value of \f$ \log_{10} \beta\gamma \f$ above which the high-energy form of the density effect may be used.
    double a;  ///< Parameter used to interpolate the density effect between the values of X0 and X1.
    double m;  ///< Parameter used to interpolate the density effect between the values of X0 and X1.
    double d0; ///< Low-energy density effect parameter, only non-zero for conducting materials.
    /// @}
};
///
/// \brief Define tdata.
///
/// Define a tdata variable for convenience.
///
///
typedef struct TDATA tdata;
///
/// \brief Structure to store range tables.
///
/// This structure contains a range table and its associated metadata.
///
struct RANGE_TABLE {
    double z1;          ///< The projectile charge.
    double a1;          ///< The projectile mass.
    short sswitch;      ///< The switch bit field.
    tdata *target;      ///< A pointer to the ::TDATA structure used in constructing the table.
    time_t timestamp;   ///< The time at which the table was created.
    double range[MAXE]; ///< The actual table of range values.
};
///
/// \brief Define range_table.
///
/// Define a range_table variable for convenience.
///
typedef struct RANGE_TABLE range_table;
///
/// \brief The range-energy table.
///
/// This external variable contains the range-energy tables.
///
range_table trange[MAXAB];
//
// Declare all the functions in crange.cpp.
//
void usage( char* executable );
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
///
/// \brief Complex logarithm of the Gamma function.
///
/// Computes the fully complex logarithm of the fully complex Gamma function.
/// Works in all portions of the complex plane, including the negative real
/// axis.
///
/// \param z A complex number.
///
/// \return \f$ \ln \Gamma(z) \f$ , a complex number.
///
/// \warning The Gamma function has poles at all integers \<= 0.
///
template<typename T> std::complex<T> complex_lngamma( const std::complex<T> &z )
{
    static T coeff[6]={76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5};
    T x, y;
    if(z.real() > 0) {
        x=z.real()-1.0;
        y=z.imag();
    } else {
        x=-z.real();
        y=-z.imag();
    }
    T r = sqrt((x+5.5)*(x+5.5)+y*y);
    T aterm1=y*log(r);
    T aterm2=(x+0.5)*atan2(y,(x+5.5))-y;
    T lterm1=(x+0.5)*log(r);
    T lterm2=-y*atan2(y,(x+5.5)) - (x+5.5) + 0.5*log(2.0*M_PI);
    T num=0.0;
    T denom=1.000000000190015;
    for(int j=1;j<7;j++){
        T fj=(T)j;
        T cterm=coeff[j-1]/((x+fj)*(x+fj)+y*y);
        num+=cterm;
        denom+=(x+fj)*cterm;
    }
    num*=-y;
    T aterm3=atan2(num,denom);
    T lterm3 = 0.5*log(num*num + denom*denom);
    std::complex<T> result(lterm1+lterm2+lterm3,aterm1+aterm2+aterm3);
    if(z.real() < 0){
        std::complex<T> lpi(log(M_PI), 0.0);
        std::complex<T> zpi(z.real()*M_PI, z.imag()*M_PI);
        result = lpi - (result + std::log(std::sin(zpi)));
    }
    return(result);
}
///
/// \brief Confluent hypergeometric function.
///
/// Computes the confluent hypergeometric function.  All input parameters
/// are complex numbers.  Uses the formula:
/// \f[ M(a,b,z) = 1 + \sum_{n=1} \frac{(a)_n}{(b)_n}\frac{z^n}{n!} , \f]
/// where
/// \f[ (x)_n \equiv \frac{\Gamma(x+n)}{\Gamma(x)} \f]
/// is the Pochhammer Symbol.
///
/// \param a First parameter of the hypergeometric function.
/// \param b Second parameter of the hypergeometric function.
/// \param z A complex number.
///
/// \return The value \f$ M(a,b,z) \f$ , a complex number.
///
/// \warning May not be stable for large values of \f$|z|\f$.
///
template<typename T> std::complex<T> complex_hyperg(const std::complex<T> &a,
                                                    const std::complex<T> &b,
                                                    const std::complex<T> &z)
{
    T dm = 0.0;
    std::complex<T> term(1.0, 0.0);
    std::complex<T> sumterm(1.0, 0.0);
    std::complex<T> previousterm;
    do {
        previousterm = term;
        dm += 1.0;
        std::complex<T> Cm(dm-1.0, 0.0);
        term = previousterm * ((a + Cm)/(b + Cm)) * (z/dm);
        sumterm += term;
    } while( std::abs(term) > 1.0e-6 && std::abs(previousterm) > 1.0e-6 );
    return(sumterm);
}
} // end namespace CRange
