///
/// \mainpage crange - The Berkeley Range-Energy Calculator
/// \author Benjamin Weaver <benjamin.weaver@nyu.edu>
/// \version 2.0.0
///
/// \copyright (C) 2001-2016 Benjamin Weaver, LGPL
///
/// \section sec-intro Introduction
///
/// Thank you for choosing the Berkeley Range-Energy Calculator.
///
/// \section sec-preinstall Pre-Installation
///
/// \subsection subsec-sysreq System Requirements
///
/// crange may require cmake to be installed, depending on how you
/// obtained the code.
///
/// crange uses Doxygen for documentation processing.
///
/// \subsection subsec-prereq Library Prerequisites
///
/// \subsubsection subsubsec-gsl GNU Scientific Library
///
/// The GNU Scientific Library (http://www.gnu.org/s/gsl/) is required.
/// This code was most recently tested with version 1.15, though any version
/// with full complex variable support should work.
///
/// \subsubsection subsubsec-iniparser iniParser
///
/// The iniParser Library (http://ndevilla.free.fr/iniparser/), version 3.0 or
/// later is strongly recommended.  The code will still compile if it is
/// not installed, but the resulting binary will not be able to read optional
/// target or switch information, so only compiled-in defaults will be
/// available.
///
/// \section sec-install Installation
///
/// See the README.md file.
///
/// \section sec-running Running
///
/// \subsection subsec-command On the Command Line
///
/// Type crange -h to get the list of command-line options.
///
/// \subsection subsec-tasks Task List
///
/// The crange program executes a list of tasks from a file.  Here's an example:
///
/// \verbinclude tasks.txt
///
/// The first line tells the calculator to compute the range (in g cm<sup>-2</sup>) of
/// uranium (Z=92,A=238) at a kinetic energy of 1200 MeV per nucleon, in an
/// aluminum target.  The second line tells the calculator to compute the
/// initial kinetic energy (in A MeV) of gold (Z=79,A=197) whose range in
/// the plastic track-etch detector CR-39 was 9.2 g cm<sup>-2</sup>.  The third line asks
/// for the final kinetic energy after passing through 9.2 g cm<sup>-2</sup> of CR-39, given
/// an initial energy of 10.6 A GeV. The fourth line tells the calculator to
/// compute dE/dx (in A MeV g<sup>-1</sup> cm<sup>2</sup>) for gold with kinetic energy 10.6 A GeV in
/// air.  The fifth line computes REL instead of dE/dx with the REL cutoff set to
/// 300 eV. The sixth line computes primary ionization.
/// Arguments which are zero (0) are dummies which are necessary for
/// place-holding.
///
/// The list of tasks may be of any length in any combination of ranges,
/// energies or dE/dx.  The order of the output will be the same as the
/// order of the input.
///
/// \subsection subsec-switch Switching Optional Effects
///
/// The switch.ini file included with the source distribution shows the
/// effects that may be turned on or off by the user.  The values in the
/// switch.ini file are the default values that are compiled into the program.
/// The user may supply a modified switch.ini file on the command line.
/// There are additional details in the file itself.
///
/// \subsection subsec-target Adding or Modifying Targets
///
/// The target.ini file included with the source distribution lists the
/// targets that are compiled in to the program by default.  Most of the target
/// data is taken from Sternheimer, Berger \& Seltzer, \cite art_rms1.
/// The definitions of the material properties are in the target.ini file.
/// The user may add additional targets or orverride existing target values
/// by supplying a different target.ini file (with the same format!) on
/// the command line.
///
/// \section sec-updates Updates
///
/// Visit http://cosmo.nyu.edu/~bw55/dedx/ .
///
/// \section sec-history History
///
///  - 2.0.0: Complete re-write in C++.
///  - 1.6.2: Minor changes relating to GitHub migration.
///  - 1.6.1: Fix version strings.
///  - 1.6.0: Make crange compatible with GNU autotools and GNU Scientific Library.
///  - 1.5.3: Simplified calculation switching and added a crange.h file
///  - 1.5.2: Increased value of MAXAB to acommodate a larger target data file
///  - 1.5.1: Fixed typo in one of the parameters of the electron capture corection
///
///
/// \file main.cpp
/// \brief Source code for main() function.
///
#include "crange.h"
///
/// \brief Main program.
///
/// Standard C/C++ main() program.
///
/// \param argc Number of command line options.
/// \param argv The command line options.
///
/// \return The exit status.
///
int main(int argc, char *argv[])
{
    bool listflag = false;
    std::string command, outputfile, targetfile, switchfile;
    int c;
    while ((c = getopt(argc, argv, ":c:hlo:s:t:V")) != -1) {
        switch (c) {
            case 'c':
                // Interpret the value of optarg as a command to the crange engine.
                command = optarg;
                break;
            case 'h':
                // Print help and exit.
                CRange::usage(argv[0]);
                return 0;
            case 'l':
                // Print the built-in target table.  Still needed?
                listflag = true;
                break;
            case 'o':
                // Print output to a file.
                outputfile = optarg;
                break;
            case 's':
                // Use filename as a switch file.
                switchfile = optarg;
                break;
            case 't':
                // Use filename as a target file.
                targetfile = optarg;
                break;
            case 'V':
                // Print version string and exit.
                CRange::version(argv[0]);
                return 0;
            case '?':
                // Unknown!
                std::cerr << "Unrecognised option: -" << (char)optopt << std::endl;
                CRange::usage(argv[0]);
                return 1;
        }
    }
    std::vector<CRange::Tdata> targets = CRange::init_target(targetfile);
    short sswitch = CRange::init_switch(switchfile);
    if (listflag) {
        for (std::vector<CRange::Tdata>::iterator it=targets.begin(); it!=targets.end(); ++it) std::cout << *it << std::endl;
        return 0;
    }
    std::vector<std::string> commands;
    if (argc-optind >= 1) {
        std::string taskfile = argv[optind];
        std::ifstream finput;
        finput.open(taskfile);
        if (!finput.is_open()) {
            std::cerr << "Error opening task file: " << taskfile << std::endl;
            return 2;
        }
        std::string line;
        while ( std::getline(finput, line) ) {
            commands.push_back(line);
        }
        finput.close();
    } else if (command.length() > 0) {
        //
        // Create a stream to hold the command.
        //
        commands.push_back(command);
    } else {
        std::cerr << "No task file specified!" << std::endl;
        return 2;
    }
    // CRange::Tdata foo = CRange::find_target("Hosta",targets);
    // std::cout << foo << std::endl;
    // std::cout << CRange::effective_charge(92.0, 950.0, foo.z2(), sswitch) << std::endl;
    for (std::vector<std::string>::iterator it=commands.begin(); it != commands.end(); ++it) {
        std::cout << *it << std::endl;
    }
    return 0;
}
