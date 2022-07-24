# SBody

[![CMake](https://github.com/Tai-Zhou/SBody/actions/workflows/cmake.yml/badge.svg)](https://github.com/Tai-Zhou/SBody/actions/workflows/cmake.yml)

The "S" in the name of SBody stands for "Small", "Some", and "Speed".

* [Library & Tool](#library--tool)
  * [CMake](#cmake)
  * [GSL](#gsl)
  * [CFITSIO](#cfitsio)
  * [indicators](#indicators)
  * [Python Packages](#python-packages)
    * [NumPy](#numpy)
* [Changelog](#changelog)
* [License](#license)

## Library & Tool
### CMake
[CMake](https://cmake.org) is an open-source, cross-platform family of tools designed to build, test and package software. SBody use CMake to get built on different platforms.
* Linux: run `apt install cmake` or download and install from [homepage](https://cmake.org).
* macOS: Install CMake via [Homebrew](https://brew.sh), with command `brew install cmake`.

### GSL
The GNU Scientific Library ([GSL](https://www.gnu.org/software/gsl/)) is a numerical library for C and C++ programmers.
* Linux: run `apt install libgsl-dev` or download and install from [homepage](https://www.gnu.org/software/gsl/).
* macOS: Install GSL via [Homebrew](https://brew.sh), with command `brew install gsl`.

### CFITSIO

[CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) is a library of C and Fortran subroutines for reading and writing data files in FITS (Flexible Image Transport System) data format. SBody use CFITSIO to produce output files in FITS.
* Linux: run `apt install libcfitsio-dev` or download and install from [homepage](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/).
* macOS: Install CFITSIO via [Homebrew](https://brew.sh), with command `brew install cfitsio`, you may need to add the library directory to `CPLUS_INCLUDE_PATH`.

### indicators
SBody uses [indicators](https://github.com/p-ranav/indicators) to show the progress.

### Python Packages
#### NumPy
SBody output can be saved in [NPY format version 2.0](https://numpy.org/devdocs/reference/generated/numpy.lib.format.html#format-version-2-0).

## Changelog
Changelog can be found [here](CHANGELOG.md).

## License
This program uses [GNU Scientific Library](https://www.gnu.org/software/gsl/) (GSL), so the program can only be distributed under the terms of the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html) (GPL).
