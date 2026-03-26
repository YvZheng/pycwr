# -*- coding: utf-8 -*-
from pathlib import Path

from setuptools import Extension, find_namespace_packages, setup

ROOT = Path(__file__).resolve().parent
PACKAGE_ROOT = ROOT / "pycwr"
RADAR_GRID_BASE = PACKAGE_ROOT / "core" / "RadarGridC"
RADAR_GRID_MODULE_BASE = Path("pycwr") / "core" / "RadarGridC"

DISTNAME = "pycwr"
AUTHOR = "pycwr developers"
AUTHOR_EMAIL = "YuZheng1206@outlook.com"
URL = "https://github.com/YvZheng/pycwr"
LICENSE = "PolyForm-Noncommercial-1.0.0"
PYTHON_REQUIRES = ">=3.9"
INSTALL_REQUIRES = [
    "numpy",
    "pandas<3",
    "pyproj",
    "scipy",
    "xarray",
    "netCDF4",
    "easydict",
]
EXTRAS_REQUIRE = {
    "full": [
        "matplotlib",
        "cartopy",
        "Flask>=3.0",
        "arm_pyart>=2.1,<3; python_version >= '3.10'",
        "xradar>=0.7; python_version >= '3.10'",
    ],
}
DESCRIPTION = "pycwr 1.0.2: high-compatibility Chinese weather radar IO, geometry, QC, plotting, and compositing"
LONG_DESCRIPTION = (ROOT / "README.md").read_text(encoding="utf-8")
PLATFORMS = ["Linux", "Mac OS-X", "Windows"]
CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Developers",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Atmospheric Science",
    "Operating System :: OS Independent",
]


def _build_extensions():
    import numpy

    include_dirs = [numpy.get_include()]
    c_source = RADAR_GRID_BASE.with_suffix(".c")
    c_module_source = RADAR_GRID_MODULE_BASE.with_suffix(".c")
    if c_source.exists():
        extension = Extension(
            "pycwr.core.RadarGridC",
            [str(c_module_source)],
            include_dirs=include_dirs,
        )
        return [extension], include_dirs

    pyx_source = RADAR_GRID_BASE.with_suffix(".pyx")
    pyx_module_source = RADAR_GRID_MODULE_BASE.with_suffix(".pyx")
    if not pyx_source.exists():
        raise RuntimeError(
            "Missing extension sources: expected pycwr/core/RadarGridC.c or "
            "pycwr/core/RadarGridC.pyx."
        )

    try:
        from Cython.Build import cythonize
    except ImportError as exc:
        raise RuntimeError(
            "Missing generated extension source: pycwr/core/RadarGridC.c. "
            "Install Cython or regenerate the C source from pycwr/core/RadarGridC.pyx "
            "before building a distribution."
        ) from exc

    extension = Extension(
        "pycwr.core.RadarGridC",
        [str(pyx_module_source)],
        include_dirs=include_dirs,
    )
    return cythonize([extension], compiler_directives={"language_level": "3"}), include_dirs


EXT_MODULES, INCLUDE_DIRS = _build_extensions()


setup(
    name=DISTNAME,
    version="1.0.2",
    author=AUTHOR,
    license=LICENSE,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    url=URL,
    platforms=PLATFORMS,
    classifiers=CLASSIFIERS,
    include_package_data=True,
    package_data={
        "pycwr": [
            "data/*",
            "draw/colormap/balance-rgb.txt",
            "GraphicalInterface/templates/*.html",
            "GraphicalInterface/static/*",
        ],
    },
    packages=find_namespace_packages(include=["pycwr", "pycwr.*"]),
    ext_modules=EXT_MODULES,
    include_dirs=INCLUDE_DIRS,
    zip_safe=False,
)
