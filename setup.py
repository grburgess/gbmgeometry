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


# # Load the package's __version__.py module as a dictionary.
# about = {}
# if not VERSION:
#     with open(os.path.join(here, NAME, "__version__.py")) as f:
#         exec(f.read(), about)
# else:
#     about["__version__"] = VERSION


# class UploadCommand(Command):
#     """Support setup.py upload."""

#     description = "Build and publish the package."
#     user_options = []

#     @staticmethod
#     def status(s):
#         """Prints things in bold."""
#         print("\033[1m{0}\033[0m".format(s))

#     def initialize_options(self):
#         pass

#     def finalize_options(self):
#         pass

#     def run(self):
#         try:
#             self.status("Removing previous builds...")
#             rmtree(os.path.join(here, "dist"))
#         except OSError:
#             pass

#         self.status("Building Source and Wheel (universal) distribution...")
#         os.system("{0} setup.py sdist bdist_wheel --universal".format(sys.executable))

#         self.status("Uploading the package to PyPI via Twine...")
#         os.system("twine upload dist/*")

#         self.status("Pushing git tags...")
#         os.system("git tag v{0}".format(about["__version__"]))
#         os.system("git push --tags")

#         sys.exit()
