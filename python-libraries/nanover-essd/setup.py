#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_namespace_packages

setup(
    name="nanover-essd",
    version="1.0.0",
    description="Extremely Simple Server Discovery, for use with NanoVer",
    author="Intangible Realities Lab",
    author_email="m.oconnor@bristol.ac.uk",
    url="https://gitlab.com/intangiblerealities/",
    packages=find_namespace_packages("src", include="nanover.*"),
    package_data={"": ["py.typed"]},
    install_requires=("netifaces2",),
    entry_points={
        "console_scripts": ["nanover-essd-list=nanover.essd.list_cli:main"],
    },
    package_dir={"": "src"},
)
