from setuptools import setup
import versioneer




setup(
    license="BSD",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    scripts = ['bin/get_trigdat', 'bin/get_posthist']
)
