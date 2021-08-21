from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import ldpc
import os

from distutils.dir_util import copy_tree
from shutil import rmtree

ldpc_path=ldpc.get_include()+'/include'

source_files=["src/bposd/bposd.pyx",
            "src/bposd/include/binary_char.c",
            "src/bposd/include/mod2sparse_extra.c",
            "src/bposd/include/sort.c",
            ldpc_path+"/mod2sparse.c"]

extension = Extension(
    name="bposd.bposd",
    sources=source_files,
    libraries=[],
    library_dirs=[],
    include_dirs=[ldpc.get_include(), ldpc_path, numpy.get_include(),"src/bposd/include"],
    extra_compile_args=['-std=c11']
    )

setup(
    python_requires='>=3.6',
    name='bposd',
    version='0.0.1',
    description='BP+OSD',
    long_description='This module provides an implementation of the belief\
        propagagation + ordered statistics decoder for quantum LDPC codes.',
    url='https://roffe.eu',
    author='Joschka Roffe',
    packages=["bposd"],
    package_dir={'':'src'},
    package_data = {'bposd': ['*.pxd']},
    include_package_data=True,
    zip_safe=False,
    ext_modules=cythonize([extension]),
    classifiers=['Development Status :: 1 - Planning'],
    install_requires=["tqdm","scipy","ldpc","numpy"]

)

