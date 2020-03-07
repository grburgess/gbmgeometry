from setuptools import setup, find_packages, Command
import os
import io
import sys
from shutil import rmtree
import versioneer

setup_requires = ["setuptools >= 30.3.0"]


setup(
    setup_requires=setup_requires,
    license="BSD",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
