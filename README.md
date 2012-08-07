vmc
===

Performs variational Monte Carlo (VMC) calcuations on lattice systems
that are of interest to the authors of this code.

This package is made up of two distinct components:

* vmc-core is the heart of the Monte Carlo calculation

* pyvmc is a higher-level package that allows setting up, controlling,
  and analysis of calculations.

Currently, the two components communicate by using JSON over a pipe.
A lot of functionality is duplicated in each component so that they
can talk to each other in the same terms.  Over time, we'd like to
remove this complexity by migrating toward Cython.

Compiling vmc-core
------------------

Requires a recent [boost](http://www.boost.org/) (headers only) and
[eigen3](http://eigen.tuxfamily.org/).  The input mechanism for
declaring calculations also requires
[jsoncpp](http://jsoncpp.sourceforge.net/) 0.6.0-rc2 or later.

A file named vmc-core/Makefile-vmc.local can be created to override
any variables in the Makefile. For example, it could say:

    EIGEN3_CFLAGS = -I/path/to/eigen3/include
    BOOST_CFLAGS = -I/path/to/boost/include
    CXX = clang++

To compile and run:

    $ cd vmc-core
    $ make && ./vmc-core < sample-input.json

It should compile on recent versions of g++, clang++, and icc.

vmc-core user documentation
---------------------------

The "vmc-core" program uses JSON for input and output.  A sample input
file is given as sample-input.json.  More complete documentation about
its structure may exist soon.

vmc-core API documentation
--------------------------

If doxygen is installed, vmc-core API documentation can be generated
by running

    $ make docs

Afterwards, HTML documentation will exist at
vmc-core/docs/generated/html/index.html

If PDF output is desired, change directory to
vmc-core/docs/generated/latex/ and run "make".

vmc-core tests
--------------

There is not standalone test framework for vmc-core, but we use
assertions liberally.  There are some tests/assertions that slow the
code significantly if enabled.  To run with these tests, compile with
the flag -DVMC_CAREFUL and run vmc-core.

If we build standalone tests for vmc-core in the future we will likely
use googletest.

Another outstanding task is to write a fuzz test for CeperleyMatrix.

pyvmc
-----

Documentation on performing a simple calculation using pyvmc is coming
soon.

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

Ways in which things are currently in flux
------------------------------------------

We are currently generalizing things so that multi-particle moves can
be made.  Instead of doing this all at once, we are making progress in
phases.  At the moment the following things are broken:

* finish_update() can only be called with a single particle move

* Renyi stuff assumes only single particle moves

* RVB wave function has not yet been updated
