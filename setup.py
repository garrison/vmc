# adapted from http://wiki.cython.org/PackageHierarchy

import sys
import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# fixme: pxd files are not automatically included as dependencies.  see
# http://www.mail-archive.com/cython-dev@codespeak.net/msg09729.html

def make_extension(ext_name, ext_libraries=()):
    return Extension(
        ext_name,
        [ext_name.replace(".", os.path.sep) + ".pyx"],
        include_dirs=(["vmc-core", "."]),
        language="c++",
        libraries=ext_libraries,
    )

extensions = [
    make_extension("pyvmc.core.lattice", ["vmc-core"]),
    make_extension("pyvmc.core.subsystem", ["vmc-core"]),
    make_extension("pyvmc.core.orbitals", ["vmc-core"]),
    make_extension("pyvmc.core.wavefunction", ["vmc-core"]),
    make_extension("pyvmc.core.measurement", ["vmc-core"]),
    make_extension("pyvmc.core.walk", ["vmc-core"]),
    make_extension("pyvmc.library.renyi", ["vmc-core"]),
    make_extension("pyvmc.core.rng", ["vmc-core"]),
    make_extension("pyvmc.core.simulation", ["vmc-core"]),
]

setup(**{
    "name": "pyvmc",
    "packages": [
        "pyvmc",
        "pyvmc.core",
    ],
    "ext_modules": extensions,
    "cmdclass": {'build_ext': build_ext},
})
