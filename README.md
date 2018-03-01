# FEniCS_Particles

Alejandro F Queiruga  
2015, 2018  
UC Berkeley, LBNL

## Intro

This is a simple library for calculating partcile-fluid interactions with FEniCS.
It was originally developed for shockwave sclae behavior with a supersonic DG solver,
but can be used for general particle-laden flows and was cleaned up for
proppant-transport problems.

## Usage and installation

This library requires C++ code which is compiled using the FEniCS' `compile_extension_module`
routine. The library will appear in your instant cache, wherever you install this repository to.

## License

Copyright (C) Alejandro Francisco Queiruga, 2015-2018

This code is released under version 3 of the GNU Lesser General Public License, as per LICENSE.txt.
