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

Note for developers: instant hashes the header file in the compiled module to determine its presence in the cache.
If you modify the cpp only, instant won't recognize the change and recompile the module.
My trick while I am developing to avoid recompiling everything else in the cache is to add random white space to the header file to change its hash value.

An example one-way coupled code is included in [examples/test.py](examples/test.py).

## License

Copyright (C) Alejandro Francisco Queiruga, 2015-2018

This code is released under version 3 of the GNU Lesser General Public License, as per LICENSE.txt.
