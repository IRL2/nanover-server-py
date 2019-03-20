#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='Narupa',
      version='1.0',
      description='Narupa python framework',
      author='Intangible Realities Lab',
      author_email='m.oconnor@bristol.ac.uk',
      url='https://gitlab.com/intangiblerealities/',
      packages=find_packages(),
      install_requires=requirements,
     )