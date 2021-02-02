#!/usr/bin/env python
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

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
      package_data={
          '': ['py.typed']
      },
      install_requires=(
            'narupa',
            'mpi4py',
            'numpy',
      ),
     )
