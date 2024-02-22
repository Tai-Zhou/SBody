# SBody

[![CMake](https://github.com/Tai-Zhou/SBody/actions/workflows/cmake.yml/badge.svg)](https://github.com/Tai-Zhou/SBody/actions/workflows/cmake.yml)

The "S" in the name of SBody stands for "Small", "Some", and "Speed".

- [Library \& Tool](#library--tool)
  - [Required](#required)
    - [CMake](#cmake)
    - [GSL](#gsl)
  - [Submodule](#submodule)
  - [Optional](#optional)
- [Changelog](#changelog)
- [License](#license)

## Library & Tool
### Required
#### CMake

[CMake](https://cmake.org) is an open-source, cross-platform family of tools designed to build, test and package software. SBody use CMake to get built on different platforms.
* Linux: run `apt install cmake` or download and install from [homepage](https://cmake.org).
* macOS: Install CMake via [Homebrew](https://brew.sh), with command `brew install cmake`.

#### GSL

The GNU Scientific Library ([GSL](https://www.gnu.org/software/gsl/)) is a numerical library for C and C++ programmers.
* Linux: run `apt install libgsl-dev` or download and install from [homepage](https://www.gnu.org/software/gsl/).
* macOS: Install GSL via [Homebrew](https://brew.sh), with command `brew install gsl`.

### Submodule

- **[indicators](https://github.com/p-ranav/indicators)**: show the progress bar.

- **[fmt](https://github.com/fmtlib/fmt)**: format the output.

### Optional

- [Doxygen](https://www.doxygen.nl): generate the documentation.
- [OpenMP](https://www.openmp.org): build parallel version.
- [Python](https://www.python.org):
  - [Breathe](https://breathe.readthedocs.io): convert the Doxygen documentation to Sphinx.
  - [NumPy](https://numpy.org): SBody output can be saved in [NPY format version 2.0](https://numpy.org/devdocs/reference/generated/numpy.lib.format.html#format-version-2-0).
  - [pybind11](https://pybind11.readthedocs.io): build Python interface.
  - [Sphinx](https://www.sphinx-doc.org): documentation shown in the webpage.

## Changelog
Changelog can be found [here](https://github.com/Tai-Zhou/SBody/blob/main/CHANGELOG.md).

## License
This program uses [GSL](https://www.gnu.org/software/gsl/), so the program can only be distributed under the terms of the [GNU General Public License](LICENSE) (GPL).
