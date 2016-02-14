# crange

## Introduction

Author: Benjamin Weaver
Email:  benjamin.weaver@nyu.edu
Homepage: http://cosmo.nyu.edu/~bw55/dedx/
Source: https://github.com/weaverba137/crange

Copyright (c) 2001-2013 by Benjamin Weaver.  This program is free software
which I release under the GNU Lesser General Public License.
Read the COPYRIGHT file for more info.  Read the INSTALL file for installation
instructions.

## Install

If you are exporting from git to create a distribution:
1. cmake .
2. make
3. make doc

If you are installing from a distribution:
1. make
2. make doc
3. make install

## Documentation

The main documentation for this package is processed by doxygen.  It is
available in PDF or HTML form.  It will be installed in
$prefix/share/doc/crange, where $prefix is /usr/local or the value
passed to the configure script with --prefix.

## TODO

Allow additional sections of the ini file to specify additional
parameters such as the size of the range table and the minimum
and maximum energy.  Sections of the ini file that do not define
switches or other configuration may be interpreted as additional
entries for the absorber table.

The pair production and bremsstrahlung computations should be updated
to reflect Sørensen's 2010 paper (see art_ahs1 in the crange.bib file).

## License

This code is released under a 3-clause BSD-style license. See the file `LICENSE.md`.
