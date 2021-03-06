///
/// \file crange.h
/// \brief Header file for crange.
///
/// This header file collects all the other header files needed to compile
/// crange, as well as all defines, function declarations, etc.
///
#ifndef _HAVE_CRANGE_H_
#define _HAVE_CRANGE_H_
//
// Standard includes here to keep crange.cpp clean
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>    // Standard (non-complex) math functions.
#include <vector>
#include <complex>  // Math functions provided here should be prefaced with std::.
#include <unistd.h> // Provides getopt(), access().
#include "config.h" // Provides version information.
#include "iniparser.h"
//
// These values define the range-energy tables
//
#define LOGTENEMIN 0.0 ///< \f$ \log_{10} E_{\mathrm{min}} \f$ Minimum energy in units of A MeV.
#define LOGTENEMAX 6.0 ///< \f$ \log_{10} E_{\mathrm{max}} \f$ Maximum energy in units of A MeV.
#define MAXE 200 ///< The number of energies in the range-energy tables.
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
///
/// \brief This type is used a lot, so good to abbreviate.
///
typedef std::complex<double> Cdouble;
///
/// \namespace CRange
/// \brief Container for most functions and classes.
///
namespace CRange
{
    ///
    /// \brief Class containing target data.
    ///
    /// This structure contains all the data related to target materials.
    ///
    class Tdata {
        private:
            static const int Ndata = 12; ///< Length of data array.
            static const std::string dnames[Ndata];  ///< Names of the data array entries.
            static const std::string comment[Ndata]; ///< Comments to be inserted when serializing.
            static const int precision[Ndata]; ///< Precision to be used when serializing.
            static const bool fixed[Ndata]; ///< Fixed or scientific notation when serializing.
            std::string _name; ///< Name of the target.
            double data[Ndata];  ///< The data.
            double _hash; ///< Computed from the name & data, used to test for equality.
            void init(void);
        public:
            Tdata();
            Tdata( const Tdata &t );
            Tdata( const std::string &n, const double d[] );
            Tdata( const char *n, const double d[] );
            Tdata( const char *n, dictionary *ini );
            void operator=( const Tdata &t );
            void print( std::ostream *o ) const;
            inline double hash (void) const { return _hash; } ///< Used in testing for equality.
            ///
            /// \name Material name
            /// @{
            ///
            const std::string& name (void) const { return _name; } ///< The target name.
            /// @}
            ///
            /// \name General parameters
            /// @{
            ///
            inline double z2  (void) const { return data[0]; } ///< The mean nuclear charge.
            inline double a2  (void) const { return data[1]; } ///< The mean atomic mass number.
            inline double iadj(void) const { return data[2]; } ///< The logarithmic mean ionization potential [eV].
            inline double rho (void) const { return data[3]; } ///< The density [g cm<sup>-3</sup>].
            inline double pla (void) const { return data[4]; } ///< The plasma frequency [eV].
            inline double etad(void) const { return data[5]; } ///< The ratio of density to density at STP for gaseous targets. Should be set to zero for non-gaseous materials.
            inline double bind(void) const { return data[6]; } ///< The total electronic binding energy [keV].
            //// @}
            ///
            /// \name Density effect parameters
            /// @{
            ///
            inline double X0  (void) const { return data[7]; } ///< Value of \f$ \log_{10} \beta\gamma \f$ at which the density effect turns on.
            inline double X1  (void) const { return data[8]; } ///< Value of \f$ \log_{10} \beta\gamma \f$ above which the high-energy form of the density effect may be used.
            inline double a   (void) const { return data[9]; } ///< Parameter used to interpolate the density effect between the values of X0 and X1.
            inline double m   (void) const { return data[10]; } ///< Parameter used to interpolate the density effect between the values of X0 and X1.
            inline double d0  (void) const { return data[11]; } ///< Low-energy density effect parameter, only non-zero for conducting materials.
            /// @}
    };
    const std::string Tdata::dnames[Tdata::Ndata] = {
        "z2", "a2", "iadj", "rho", "pla", "etad",
        "bind", "X0", "X1", "a", "m", "d0"};
    const std::string Tdata::comment[Tdata::Ndata] = {
        "Mean nuclear charge",
        "Mean nuclear mass",
        "Logarithmic mean ionization potential [eV]",
        "Density [g cm^-3]",
        "Plasma frequency [eV]",
        "Ratio of density to density at STP for gasses (zero for everything else)",
        "Total electronic binding energy [eV]",
        "Density effect turn-on value",
        "Density effect asymptotic bound",
        "Density effect interpolation parameter",
        "Density effect interpolation parameter",
        "Low energy density effect value (zero for everything but conductors)"};
    const int Tdata::precision[Tdata::Ndata] = {
        3, 3, 1, 4, 4, 1, 4, 4, 4, 5, 4, 2};
    const bool Tdata::fixed[Tdata::Ndata] = {
        true, true, true, false, false, true,
        false, true, true, true, true, true};
    ///
    /// \brief Class to store range tables.
    ///
    /// This class contains a range table and its associated metadata.
    ///
    class RangeTable {
        private:
            double range[MAXE]; ///< The actual table of range values.
        public:
            double z1;          ///< The projectile charge.
            double a1;          ///< The projectile mass.
            short sswitch;      ///< The switch bit field.
            Tdata target;       ///< The Tdata class used in constructing the table.
            RangeTable();
            RangeTable(const RangeTable &t);
            void operator=(const RangeTable &t);
            RangeTable(double z, double a, short s, Tdata &t);
            double interpolate_range( double e );
            double interpolate_energy( double e, double r0 );
    };
    //
    // These functions compute physical effects related to dE/dx.
    //
    double effective_charge( double z0, double e1, double z2, short sswitch );
    double djdx( double e1, double z0, double I0, double f0, double K, short sswitch, Tdata &target);
    double dedx( double e1, double rel0, double z0, double a1, short sswitch, Tdata &target );
    double delta( double g, Tdata &target );
    double olddelta( double g, Tdata &target );
    double bma( double z1, double b );
    double relbloch( double z12, double b1, double lambda, double theta0 );
    double lindhard( double zz, double aa, double bb, short sswitch );
    double Fbrems( double x );
    //
    // These functions integrate dE/dx to get range or energy.
    //
    double integrate_dedx( int i, double z, double a, short s, Tdata &t);
    double range( double e, double z1, double a1, short sswitch, Tdata &target, std::vector<RangeTable> &rt );
    double qrange( double e, double z1, double a1, short sswitch, Tdata &target );
    double benton( double e, double z1, double a1, Tdata &target );
    double renergy( double e, double r0, double z1, double a1, short sswitch, Tdata &target, std::vector<RangeTable> &rt );
    //
    // Utility funtions for IO, command processing, etc.
    //
    void usage( char* executable );
    void version( char* executable );
    std::vector<std::string> process( std::vector<std::string> &commands, short sswitch, std::vector<Tdata> &targets );
    short init_switch( const std::string &switchfile );
    short init_switch( const char *switchfile );
    std::vector<Tdata> init_target( const char *targetfile );
    std::vector<Tdata> init_target( const std::string &targetfile );
    Tdata find_target(const char *name, std::vector<Tdata> &targets);
    Tdata find_target(const std::string &name, std::vector<Tdata> &targets);
    std::vector<Tdata>& default_target(std::vector<Tdata> &extratargets);
    ///
    /// \brief Returns the energy corresponding to a value in a range table.
    ///
    /// This utility returns an energy value from a (virtual) vector containing
    /// A logarithmically uniform distribution of energies between a minimum
    /// and maximum energy (defined by #LOGTENEMIN and #LOGTENEMAX),
    /// with a number of entries given by #MAXE.
    ///
    /// \param i The index of the vector.
    ///
    /// \return The \em i th energy in A MeV.
    ///
    inline double energy_table( int i ) { return exp(M_LN10*(LOGTENEMIN + ((double)i)*(LOGTENEMAX-LOGTENEMIN)/(MAXE - 1.0))); }
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
    template<typename T> std::complex<T> lngamma( const std::complex<T> &z )
    {
        const static T coeff[6] = {76.18009172947146,
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
            result = lpi - (result + std::log(std::sin(M_PI*z)));
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
    template<typename T> std::complex<T> hyperg(const std::complex<T> &a,
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
//
// Requires declaration outside the namespace!
//
///
/// \brief Allow CRange::Tdata objects in output streams.
///
/// This allows a CRange::Tdata object to write itself into an output stream.
///
/// \param o The output stream (left-hand side).
/// \param t The CRange::Tdata object (right-hand side).
///
/// \return A reference to the output stream.
///
std::ostream &operator<< (std::ostream &o, const CRange::Tdata &t) { t.print(&o); return o; }
///
/// \brief Equality of two CRange::Tdata objects.
///
/// Allow two CRange::Tdata objects to be tested for equality.
///
/// \param lhs Left-hand side.
/// \param rhs Right-hand side.
///
/// \return Are the objects equal?
///
bool operator== (const CRange::Tdata &lhs, const CRange::Tdata &rhs) { return lhs.hash() == rhs.hash(); }
#endif // end ifndef _HAVE_CRANGE_H_
