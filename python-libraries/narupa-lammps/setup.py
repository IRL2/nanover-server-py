#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_namespace_packages

setup(name='narupa-lammps',
      version='0.1.0',
      description='LAMMPS integration for Narupa',
      author='Intangible Realities Lab',
      author_email='simonbennie@gmail.com',
      url='https://gitlab.com/intangiblerealities/',
      packages=find_namespace_packages('src', include='narupa.*'),
      package_dir={'': 'src'},
      install_requires=(
            'narupa',
            'mpi4py',
            'numpy',
      ),
     )
