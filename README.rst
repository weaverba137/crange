******
CRange
******

Introduction
------------

:Author: Benjamin Weaver
:Email: baweaver@lbl.gov
:Homepage: https://weaverba137.github.io/crange
:Source: https://github.com/weaverba137/crange

Install
-------

If you're installing from the tarball (provided by GitHub or otherwise):

1. ``cmake .``
2. ``make install``

This will install everything in the ``/usr/local`` path or the equivalent,
depending on your system.  To customize the install path, just do, *e.g.*:

1. ``cmake -DCMAKE_INSTALL_PREFIX=$HOME/CRange/2.0.0 .``
2. ``make install``

Depending on your install path, you may need to set the environment variables
``PATH`` and ``LD_LIBRARY_PATH`` (``DYLD_LIBRARY_PATH``) to point to the
executable and the shared library, respectively.

Documentation
-------------

The main documentation for this package is processed by `Doxygen`_.  It is
available in HTML form.  It will be installed in
``$prefix/share/doc/CRange/html``, where ``$prefix`` is ``/usr/local`` or the value
passed to the cmake command with ``-DCMAKE_INSTALL_PREFIX``.

.. _`Doxygen`: http://www.doxygen.org

TODO
----

* Allow additional sections of the ini file to specify additional
  parameters such as the size of the range table and the minimum
  and maximum energy.  Sections of the ini file that do not define
  switches or other configuration may be interpreted as additional
  entries for the absorber table.
* The pair production and bremsstrahlung computations should be updated
  to reflect Sørensen's 2010 paper (see art_ahs1 in the crange.bib file).

License
-------

This code is released under a 3-clause BSD-style license. See the file ``LICENSE.rst``.

This package includes code from the `iniparser library`_, version 4.0,
Copyright (c) 2000-2011 by Nicolas Devillard.  See the file ``INIPARSER_LICENSE.rst``.

.. _`iniparser library`: https://github.com/ndevilla/iniparser
