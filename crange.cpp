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
/// \brief Print usage message.
///
/// Prints a usage message on STDERR.
///
/// \param executable Name of the program, usually argv[0].
///
void CRange::usage( char *executable )
{
    std::cerr << "usage: " << executable << " [-c COMMAND] [-h] [-l] [-o FILE] [-s switch.ini] [-t target.ini] <task file>" << std::endl;
    std::cerr << "       -c COMMAND    = Execute this one-line command instead of reading it from a file." << std::endl;
    std::cerr << "       -h            = Print this help message and exit." << std::endl;
    std::cerr << "       -l            = Print the built-in target table and exit." << std::endl;
    std::cerr << "       -o FILE       = Write to this file instead of standard output." << std::endl;
    std::cerr << "       -s switch.ini = Override the default switch values by reading this file." << std::endl;
    std::cerr << "       -t target.ini = Override the default target values by reading this file." << std::endl;
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
