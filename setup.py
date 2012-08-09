# adapted from http://wiki.cython.org/PackageHierarchy

import sys
import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

def make_extension(ext_name, ext_extras=(), ext_libraries=()):
    ext_path = ext_name.replace(".", os.path.sep) + ".pyx"

    return Extension(
        ext_name,
        ([ext_path] + [os.path.join("vmc-core", e) for e in ext_extras]),
        include_dirs=(["vmc-core", "."]),
        language="c++",
        #define_macros=[],
        extra_compile_args=["-O3", "-Wall", "-Wno-self-assign"],
        #extra_link_args = ['-g'],
        libraries=ext_libraries,
    )

extensions = [
    make_extension("pyvmc.core.lattice", ["Lattice.cpp"]),
]

setup(**{
    "name": "pyvmc",
    "packages": [
        "pyvmc", "pyvmc.core"
    ],
    "ext_modules": extensions,
    "cmdclass": {'build_ext': build_ext},
})