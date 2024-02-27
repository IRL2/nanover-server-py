#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_namespace_packages

setup(
    name="nanover-lammps",
    version="0.1.0",
    description="LAMMPS integration for NanoVer",
    author="Intangible Realities Lab",
    author_email="simonbennie@gmail.com",
    url="https://gitlab.com/intangiblerealities/",
    packages=find_namespace_packages("src", include="nanover.*"),
    package_dir={"": "src"},
    package_data={"": ["py.typed"]},
    install_requires=(
        "nanover",
        "mpi4py",
        "numpy",
    ),
)
