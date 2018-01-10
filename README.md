# FEM
Generalized, arbitrary order continuous finite element code on triangles

# Building
The only dependency is [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) (either the serial or multi threaded openMP versions will work). This can be avoided by not using the LU solver and instead using Eigen's Cholesky solver. Currently, there isn't an alternative for directly solving unsymmetric matrices (such as for Navier Stokes). 

To build, copy make.inc.in to make.inc and edit the settings. Only one of SUPERLU or PSUPERLU can be defined (can use serial or multithreaded SuperLU as there is an issue with redefining blas). From there, executables can be built with 
```bash 
  make name_of_program.exe 
```
where name_of_program.cpp is the source file with int main() in it. All programs in the `exec` directory take a mesh name as an argument. For example, running `stokes.exe` on the cavity mesh: 
```bash
./stokes.exe cavity
``` 

# Meshing 
`gmsh` was used to create all the meshes in the `mesh/` directory. The python program `mesh_converter.py` converts the `gmsh` ASCII format to a simpler format this code can read in. Materials are specified by setting physical surfaces in `gmsh` and boundary conditions can be inputted by setting physical lines to the correct numbers. The way boundary conditions are handled is in the constructor of the `FEGrid` class. 
