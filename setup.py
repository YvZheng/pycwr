# -*- coding: utf-8 -*-
"""
suggested by bugsuse(https://github.com/bugsuse)
"""

from setuptools import find_packages, setup
import sys, os
from Cython.Build import cythonize
import numpy

parent_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(parent_dir)

DISTNAME = "pycwr"
AUTHOR = "pycwr developers"
AUTHOR_EMAIL = "YuZheng1206@outlook.com"
URL = "https://github.com/YvZheng/pycwr"
LICENSE='MIT'
PYTHON_REQUIRES = ">=3.6"
INSTALL_REQUIRES = ["matplotlib", "cython", "pyproj", "Cartopy", "xarray","numpy",\
                    "scipy", "pandas", "PyQt5", "netCDF4", 'easydict']
DESCRIPTION = "China Weather Radar tools"
LONG_DESCRIPTION = """The Weather Radar Toolkit, support most of China's radar formats
(WSR98D, CINRAD/SA/SB/CB, CINRAD/CC/CCJ, CINRAD/SC/CD)"""
PLATFORMS = ["Linux", "Mac OS-X", "Windows"]
CLASSIFIERS = [
    'Development Status :: 1 - Planning',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Operating System :: POSIX :: Linux',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']
setup(
    name=DISTNAME,
    version="0.3.1",
    author=AUTHOR,
    license=LICENSE,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    url=URL,
    platforms=PLATFORMS,
    classifiers=CLASSIFIERS,
    include_package_data = True,
    packages=find_packages(parent_dir),
    ext_modules=cythonize('./pycwr/core/RadarGridC.pyx'),
    include_dirs=[numpy.get_include()]
)

