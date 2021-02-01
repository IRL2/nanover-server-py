#!/usr/bin/env python
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from distutils.core import setup
from setuptools import find_namespace_packages

setup(name='narupa-essd',
      version='1.0.0',
      description='Extremely Simple Server Discovery, for use with Narupa',
      author='Intangible Realities Lab',
      author_email='m.oconnor@bristol.ac.uk',
      url='https://gitlab.com/intangiblerealities/',
      packages=find_namespace_packages('src', include='narupa.*'),
      package_data={
          '': ['py.typed']
      },
      install_requires=(
          'netifaces',
      ),
      entry_points={
          'console_scripts': ['narupa-essd-list=narupa.essd.list_cli:main'],
      },
      package_dir={'': 'src'},
)
