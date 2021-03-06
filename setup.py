# adapted from http://wiki.cython.org/PackageHierarchy

import sys
import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# Use parallel compilation on this number of cores.
nthreads = int(os.environ.get('COMPILE_NTHREADS', 0))

def make_extension(ext_name, ext_libraries=(), is_directory=False):
    ext_path = ext_name
    if is_directory:
        ext_path += ".__init__"
    return Extension(
        ext_name,
        [ext_path.replace(".", os.path.sep) + ".pyx"],
        include_dirs=(["vmc-core", "."]),
        language="c++",
        libraries=ext_libraries,
        extra_compile_args=['-std=c++11'],
        #depends=["vmc-core/libvmc-core.so"],
    )

extensions = [
    make_extension("pyvmc.core.lattice", ["vmc-core"]),
    make_extension("pyvmc.core.boundary_conditions", ["vmc-core"]),
    make_extension("pyvmc.core.subsystem", ["vmc-core"]),
    make_extension("pyvmc.core.orbitals", ["vmc-core"]),
    make_extension("pyvmc.core.wavefunction", ["vmc-core"]),
    make_extension("pyvmc.core.estimate", ["vmc-core"]),
    make_extension("pyvmc.core.measurement", ["vmc-core"]),
    make_extension("pyvmc.core.walk", ["vmc-core"]),
    make_extension("pyvmc.core.rng", ["vmc-core"]),
    make_extension("pyvmc.core.simulation", ["vmc-core"]),
    make_extension("pyvmc.measurements", ["vmc-core"], is_directory=True),
    make_extension("pyvmc.library.renyi", ["vmc-core"]),
    make_extension("pyvmc.library.dbl", ["vmc-core"]),
    make_extension("pyvmc.library.dmetal", ["vmc-core"]),
    make_extension("pyvmc.library.bcs", ["vmc-core"]),
]

setup(**{
    "name": "pyvmc",
    "packages": [
        "pyvmc",
        "pyvmc.core",
    ],
    "ext_modules": cythonize(extensions, nthreads=nthreads),
    "cmdclass": {'build_ext': build_ext},
})
