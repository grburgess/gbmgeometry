from setuptools import setup


setup(

    name="gbmgeometry",
    packages=['gbmgeometry'],
    version='v0.1',
    description=' Geometry calculations for Fermi GBM ',
    author='J. Michael Burgess',
    author_email='jmichaelburgess@gmail.com',

    requires=[
        'numpy',
        'matplotlib'
        'astropy'
    ],

    
)
