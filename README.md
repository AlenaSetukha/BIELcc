# BIELCC

Library for E/D field calculation by boundary integral equetion method. 

## Description
The library implements functionality for calculating electric and magnetic fields using the boundary integral equation method.

The library allows solving problems of electromagnetic wave scattering on perfectly conducting surfaces, homogeneous dielectrics, and plasmonic particles using this method.


## Library structure
- electrodynamics
    - the basic integral operators used in the Stratton-Chu formulas for expressing the electric and magnetic fields in space: K, K0, R
    - library of basic integrands for the potentials of the simple and double layers of the Laplace and Helmoglitz equations: kernel_lib
    - classes of materials and field radiation sources: material, radiation source
    - numerical parameters of the problem and the space of fundamental electromagnetic constants: num_param, constants
    - calculation of surface currents based on local grid cell bases (piecewise constant approximation and collocation method): currents
    - calculation of the field around an ideal conductor: get_field
- geometry
    - basic operations of analytic geometry, determination of the position of a point, geometric constants: element_geom, get area, tolerance_constants

- integrators
    - numerical integrators over cells (quadrilateral and triangular): integral_universal, integral_universal_seg
    - basic analytical integrals: integrals_analytic
    - parameters of integrators and integrands: integral_param, kernel_param

- io
    - functions for reading and writing data in the format vtk: io_vtk
    - functions for reading and writing data in the format dat: io_aereco
- mesh
    - classes for constructing geometry: point, segment, cell
    - surface mesh class: surf_mesh
    - surface mesh class close to format .vtk: rwg_surf_mesh
- examples:
    - mesh: examples of creating / using surface meshes from files .vtk, .dat
    - integrators: examples of calling functions of basic numerical integrators
    - ED_Operators: examples of calling functions of the basic integral operators of the Stratton-Chu formulas
    - field: example of calculating the electrodynamic field around a perfectly conducting surface using the Stratton-Chu formulas


## How to install library
``` bash
git clone https://github.com/AlenaSetukha/BIELcc.git
cd bielcc
mkdir build &amp;&amp; cd build
cmake ..
make