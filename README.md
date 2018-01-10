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

# Documentation 
Documentation can be built by running 
```bash
make docs
``` 
which will run doxygen on the doxygen setup file in the `docs/` directory. 

# Meshing 
`gmsh` was used to create all the meshes in the `mesh/` directory. The python program `mesh_converter.py` converts the `gmsh` ASCII format to a simpler format this code can read in. Materials are specified by setting physical surfaces in `gmsh` and boundary conditions can be inputted by setting physical lines to the correct numbers. The way boundary conditions are handled is in the constructor of the `FEGrid` class. 

# Source Files
## fem 
Contains all the souce files related to FEM. `FEGrid` stores a collection of `Elements` and `Nodes` as specified by the mesh. All elemental integrals are defined in the `Elements` class which stores a `Basis` objects to evaluate the basis functions and the affine transformation stuff (Jacobian, inverse Jacobian, gradients of basis functions, etc). The basis functions themselves are `Poly2D` objects which simply evaluate a 2D complete polynomial. Gauss Quadrature routines are also provided here. 

## field
Files for the `Fields` class. This class is responsible for mapping unknowns to their global index in the left hand side matrix and right hand side vector. This class is useful for problems where multiple variables are concatenated into one system such as in Navier Stokes equations where x and y velocity and pressure are simulataneously solved for. An example use of this class is:
```c++
int global_id = fields[a_node_number]["a_variable_name"]; 
```
In this way, the typical FEM assembly process can be done. The fields class converts a global FEM node number and a variable name to the corresponding matrix entry. 

## materials
Files related to storing material information. Creating materials is as simple as:
```c++
Materials mat; 
mat("material_1", "variable_name") = 5; 
mat("material_2", "variable_name") = 1; 
``` 
This information can be extracted with: 
```c++
double value = mat("material_1", "variable_name"); 
``` 

## operators
Left and right hand side assembly classes. Takes in an FEGrid and Materials object and builds the system to be solved. 

## solvers 
A slew of linear system solvers based around a custom `SparseMatrix` compressed column sparse matrix representation class. External solvers are called by converting a `SparseMatrix` to the format required by Eigen and SuperLU. 
