#!/usr/bin/env python


from distutils.core import setup
from setuptools import find_namespace_packages

setup(
    name="nanover-mdanalysis",
    version="0.1.0",
    description="MDAnalysis integration for NanoVer",
    author="Intangible Realities Lab",
    author_email="m.oconnor@bristol.ac.uk",
    url="https://gitlab.com/intangiblerealities/",
    packages=find_namespace_packages("src", include="nanover.*"),
    package_dir={"": "src"},
    package_data={"": ["py.typed"]},
    install_requires=(
        "nanover",
        "MDAnalysis",
    ),
)
