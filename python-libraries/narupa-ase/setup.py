#!/usr/bin/env python
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from distutils.core import setup
from setuptools import find_namespace_packages

setup(name='narupa-ase',
      version='0.1.0',
      description='ASE integration for Narupa',
      author='Intangible Realities Lab',
      author_email='m.oconnor@bristol.ac.uk',
      url='https://gitlab.com/intangiblerealities/',
      packages=find_namespace_packages('src', include='narupa.*'),
      package_dir={'': 'src'},
      package_data={
          '': ['py.typed']
      },
      install_requires=(
          'narupa',
          'ase>=3.21',
      ),
      entry_points={
          'console_scripts': ['narupa-omm-ase=narupa.ase.openmm.cli:main'],
      },
      )
