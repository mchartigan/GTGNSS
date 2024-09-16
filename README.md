# GTGNSS

![Static Badge](https://img.shields.io/badge/language-MATLAB-orange)

A MATLAB library for orbital dynamics, propagation, and filtering in the context of GNSS navigation.

## Requirements

1. This library requires installation of NASA JPL's SPICE utility, published through the Navigation and Ancillary Information Facility (NAIF). The MATLAB version can be installed through [the JPL NAIF website](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html).
2. Many scripts use SPICE to access data and ephemerides of celestial bodies, which stores and accesses information in [kernels](https://naif.jpl.nasa.gov/naif/spiceconcept.html). Kernels for both generic solar system data and past flight projects are available [through the NAIF website](https://naif.jpl.nasa.gov/naif/data.html) as well.

## Installation

Once the repository is cloned, add the entire `GTGNSS` directory and its subdirectories to the path. This may be done in MATLAB with the command:

```matlab
addpath(genpath('path/to/repository'))
```

## Use

Each subdirectory has modules and/or functions that can assist with orbit propagation or navigation filtering. Each module and function has a header describing the inputs (and outputs, if unclear) required to use each utility. Specific data sources are provided in the `res` directory when applicable, such as spherical harmonic coefficients for the Earth and Moon. For help information, use the command `help function-or-module-name`. If more specificity is needed, use `edit function-or-module-name` and peruse the in-line code documentation.

Please notify @mchartigan of any deficiencies in documentation or broken functionality, as this repository is a port from previous codebases.
