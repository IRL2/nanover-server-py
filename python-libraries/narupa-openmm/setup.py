#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_namespace_packages

setup(name='narupa-openmm',
      version='0.1.0',
      description='OpenMM server for Narupa',
      author='Intangible Realities Lab',
      author_email='m.oconnor@bristol.ac.uk',
      url='https://gitlab.com/intangiblerealities/',
      packages=find_namespace_packages(include='narupa.*'),
      requires=(
            'narupa',
            'openmm',
      ),
     )