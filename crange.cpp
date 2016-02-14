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
/// crange should run on any POSIX-based system that has the appropriate
/// libraries installed (see \ref subsec-prereq "Prerequisites").
///
/// \subsection subsec-prereq Library Prerequisites
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
/// See the INSTALL file.
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
///
/// \file crange.cpp
/// \brief Source code for crange.
///
/// This file contains all source code for the crange executable.
///
//
// Include headers for variable and function declarations, &c..
//
#include <crange.h>
///
/// \brief Main crange program.
///
/// Standard C/C++ main() program.
///
/// \param argc Number of command line options.
/// \param argv The command line options.
///
/// \return The exit status.
///
int main( int argc, char **argv )
{
    FILE *finput,*foutput;
    tdata *extratargets, *listdummy;
    short sswitch;
    int have_switch=0, have_target=0, have_command=0, have_output=0;
    char inputname[50];
    char *switchfile, *targetfile, *command, *outputname;
    int errflag=0, listflag=0, fd=-1;
    char tempfilename[15] = "";
    static char list[5] = "List";
    extern int errno; /* From errno.h */
    extern char *optarg; /* External variable used by getopt(). */
    extern int optind, optopt; /* External variable used by getopt(). */
    /* End declarations */
    while((int c=getopt(argc, argv, ":c:hlo:s:t:")) != -1) {
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
            std::cerr << "Option -" << (char)optopt << " requires an operand." << std::endl;
            errflag++;
            break;
        case '?':
            std::cerr << "Unrecognised option: -" << (char)optopt << std::endl;
            errflag++;
        }
    }
    if (errflag) {
        CRange::usage(argv[0])
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
            std::cerr << "Error opening task file!" << std::endl;
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
        std::cerr << "No task file specified!" << std::endl;
        return(2);
    }
    if(have_output) {
        foutput=fopen(outputname, "w");
        if (foutput==NULL) {
            std::cerr << "Error opening output file!" << std::endl;
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
