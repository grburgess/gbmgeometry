from setuptools import setup, find_packages

setup(

    name="gbmgeometry",
    packages=find_packages(),
    version='v0.1',
    description=' Geometry calculations for Fermi GBM ',
    author='J. Michael Burgess',
    author_email='jmichaelburgess@gmail.com',
    requires=[
        'numpy',
        'matplotlib',
        'astropy'
    ]


)
