#!/usr/bin/env python


from distutils.core import setup
from setuptools import find_namespace_packages

setup(
    name="nanover-jupyter",
    version="0.1.0",
    description="NGLView python client for NanoVer",
    author="Intangible Realities Lab",
    author_email="harry.stroud@usc.es",
    url="https://github.com/IRL2/",
    packages=find_namespace_packages("src", include="nanover.*"),
    package_dir={"": "src"},
    package_data={"": ["py.typed"]},
    install_requires=(
        "nanover",
        "MDAnalysis",
        "numpy",
        "nglview",
    ),
)
