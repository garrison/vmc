vmc
===

Performs variational Monte Carlo (VMC) calcuations on lattice systems
that are of interest to the authors of this code.

This package is made up of two distinct components:

* vmc-core is the heart of the Monte Carlo calculation, and is written
  in C++ for speed

* pyvmc is a higher-level package that allows setting up, controlling,
  and analysis of calculations.  It is written in python for
  ease-of-programming/use and flexibility, and communicates with the
  C++ code by using [Cython](http://cython.org/).

Compiling vmc-core
------------------

[premake4](http://industriousone.com/premake) is required for the
build.

Requires a recent [boost](http://www.boost.org/) (headers only) and
[eigen3](http://eigen.tuxfamily.org/).

To compile and run (using e.g. clang++):

    $ cd vmc-core
    $ premake4 --os=linux gmake && CXX=clang++ INCLUDES="-I/path/to/eigen3/include -I/path/to/boost/include" make

Run "premake4 --help" for a list of options.  Note that it is possible
for premake to hook into a different build system instead of using
"make."

The code should compile on recent versions of g++, clang++, and icc.

vmc-core API documentation
--------------------------

If [doxygen](http://www.doxygen.org/) is installed, vmc-core API
documentation can be generated by running

    $ cd vmc-core/
    $ mkdir -p docs/
    $ doxygen

Afterwards, HTML documentation will exist at
vmc-core/docs/generated/html/index.html

If PDF output is desired, change directory to
vmc-core/docs/generated/latex/ and run "make".

vmc-core tests
--------------

There is no standalone test framework for vmc-core, but we use
assertions liberally.  There are some tests/assertions that slow the
code significantly if enabled.  To run with these tests, compile with
the flag -DVMC_CAREFUL by running:

    $ make config=careful

Be sure to pass all the appropriate environment variables, given
above, to "make."

If we build standalone tests for vmc-core in the future we will likely
use googletest.

Another outstanding task is to write a fuzz test for CeperleyMatrix.

pyvmc
-----

Documentation on performing a simple calculation using pyvmc is coming
soon.

Compiling pyvmc
---------------

First make sure a recent python2 is installed, along with virtualenv.
The following uses clang++ for both compiling (CC) and linking (CXX),
and assumes certain (typical) directories for boost and eigen.

The "pip" command below compiles all the python dependencies given in
requirements.txt.  Some of these have their own dependencies, which
must be installed.

    $ virtualenv venv
    $ source venv/bin/activate
    $ pip install -r requirements.txt
    $ CC=clang++ CXX=clang++ python setup.py build_ext -I/usr/include/boost:/usr/include/eigen3 -Lvmc-core --inplace

Sometimes you must pass the option "--force" to the build_ext command,
since [dependencies are not tracked properly when pxd files are
modified](http://www.mail-archive.com/cython-dev@codespeak.net/msg09729.html).

LD_LIBRARY_PATH
---------------

After having followed the instructions above, things will have
compiled but the operating system will not know where to look for the
vmc-core library when python needs to load it.  It can either be
copied or symlinked to e.g. /usr/local/lib, or the environment
variables can be set so that the library path is searched.  For
instance, in the bash shell on Linux:

    $ export LD_LIBRARY_PATH=/path/to/vmc/vmc-core

It is recommended to construct two scripts, one for building and one
for executing python, to cut down on the effort needed to build and
run things.  Otherwise, the environment variables and commands can get
out of hand.

pyvmc tests
-----------

py.test is used to perform unit tests of pyvmc.

    $ PYTHONPATH=. py.test pyvmc

There are also Monte Carlo tests that use both pyvmc and vmc-core to
perform calculations and compare to previous (or known) quantities.
These are located in pyvmc/mc_tests/ and can be executed by running:

    $ python pyvmc/mc_tests/__init__.py

Caveats
-------

Only wavefunctions whose (slave) particles obey Pauli exclusion are
currently supported (i.e. fermions and hard-core bosons).

Only wavefunctions with a definite particle number are supported.

Also, only wavefunctions that have a definite $S_z$ are supported.
Otherwise, it is not entirely accurate to think of spin-up and
spin-down particles as separate species (as the code does), since
annihilation operators anticommute even if the spins of the two
operators are different.

Another way of phrasing the above two requirements is that moves and
operators must consist only of SiteHop's.

Ways in which things are currently in flux
------------------------------------------

At the moment the following things are broken:

* Renyi stuff uses only single particle moves, even on wavefunctions
  where that doesn't work well.

* BCS wavefunction has not been tested.

* Projected Fermi sea does not yet use multi-particle moves so does
  not work at half filling.

* Non-Bravais lattices have never been tested, and may not be fully
  implemented in the python layer.

* Cylindrical boundary conditions are not supported.  Also, boundary
  conditions in python are the reciprocal of their representation in
  C++, which should be fixed at some point.
