# A C++ Code for Solving a Contact Problem with a Slipping Fault
## S. Galati
### A.A. 2024-2025

The project aims to solve a contact mechanics problem using different methods for imposing both the non-penetrability and the Coulomb friction condition at the interface representing the fault. The implemented code allows to solve 3-dimensional problems where the geometry is a parallelepiped bulk which is entirely cut by a fault intersecting the top and the bottom boundary, possibly with a small inclination. The codebase is compiled with \texttt{CMake} and is used to generate a small static library \texttt{mycontactlib.a} that is used by the {main.cpp} source code.
The FEM framework is given by GetFEM++, a very powerful library for solving many FEM problems. In order to compile and execute the code, the user will need the following libraries installed:
- getfem
- mumps (sequential version): for the linear algebra procedures

The implemented code allows to solve 3-dimensional problems where the geometry is a parallelepiped bulk which is entirely cut by a fault intersecting the top and the bottom boundary, possibly with a small inclination. The codebase is compiled with `CMake` and is used to generate a small static library `mycontactlib.a` that is used by the `main.cpp` source code.
It supports Gmsh integration for the generation of inclined and/or unstructured meshes. The main library used is GetFEM++ with an internal binding to some common linear algebra libraries. Before compiling the code, the user needs to have those libraries installed in their system directories.

# Notes for compiling the code
To install GetFEM++, first download the last stable version and unpack it:
```bash
    wget http://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-5.4.4.tar.gz
    tar xzf getfem-5.4.4.tar.gz
    cd getfem-5.4.4/
```

Before starting the configuration process, make sure that the following libraries are installed on your system:
  - `Qhull`: used for Delaunay triangulations;
  - `SuperLU`: A direct solver for large sparse linear systems;
  - `Mumps`: A direct solver for large sparse linear systems. Be sure to have the sequential version installed;
  - `Lapack`, `BLAS`.


Then configure GetFEM++ with:
```bash
    ./configure --enable-shared --enable-mumps --enable-superlu
```

This supposes that the used libraries are installed in system directories. If you want to use different libraries installed in non-standard paths, you can specify it with the options:
  - `--with-blas=<lib>`             use BLAS library <lib>
  - `--with-superlu=<lib>`          use SuperLU library <lib>
  - `--with-superlu-include-dir`    directory in which the
                                    superlu/sl*.h or just sl*.h   headers can be found
  - `--with-mumps=<lib>`            use MUMPS library <lib>
  - `--with-mumps-include-dir`


For more information on the configuration, run
```bash
    ./configure --help
```

To build `getfem`:
```bash
    gmake -j<nprocs>
```
You can optionally check that the installation was successful:
```bash
    gmake check
```
Finally, to install the library run:
```bash
    sudo gmake install
```

Since the code supports mesh generation via the open-source software Gmsh, you may need to install that if you plan to use it. Other optional packages are `doxygen` and `graphviz` for generating documentation.

Once you have installed the GetFEM++ library, go to the \texttt{PROJECT\_ROOT} and run:
```bash
    mkdir build && cd build
    cmake ..
```
with possible options (all defaulted to OFF):
  - `-DEXAMPLES_VERBOSE=ON`: to be verbose when testing the examples 
  - `-DEXAMPLES_USE_GMSH=ON`: to test the examples using a mesh generated via gmsh

and then:
```bash
    make -j<nprocs>
```
Note: If the `build` directory gets automatically populated by CMake due to your IDE settings when you restart the IDE, remove the directory and create a new one with the given command

Optionally, you can execute the example both to check that the installation of the library was successful and to check that the output is correct, by doing:
```bash
    ctest
```
To see all the possible options run `ctest --help`.

To clean the result of compilation, you can do
```bash
    make clean-all
```
Finally, to generate the documentation, run
```bash
    make doc    
```


### Source code
The main implemented classes are:
`Params`:
A struct holding all the parameters contained in `data.pot`.
- It implements a function for reading the parameters from a GetPot datafile, including some options passed at runtime.
- Overloads the streaming operator.

`Mesh` and `MeshBuilderStrategy`: the first is just a wrapper for the `getfem::mesh` object, which is built using the Builder in `MeshBuilder.hpp`. It allows direct access to the underlying regions definied within the mesh (Bulks, Fault and Boundary regions) and encapsulates the logic for built-in construction or Gmsh import.

`FEMManager`: to manage the `getfem::mesh_fem` objects.

`BCHandler` and `BC`: to handle the boundary conditions. The supported BC types are Dirichlet, Neumann and normal Dirichlet. In the data.pot:
- region< type > = [1 4 5]: Dirichlet is imposed on boundary regions with tags 1,4 and 5
- bd< type >< i > = {x[1]+sqrt(x[0]),x[2],3.5*t^2}: the < type > to be imposed on the i-th region given in region< type >

`ContactProblem` and `ContactEnforcementStrategy`: the main class handling the problem. It implements methods to:
- initialize the problem: set BCs, set FEM, set integration method
- prepare the assemble: use GFWL of getfem to instruct the solver on how to assemble the matrices at each iteration, choosing the method through a strategy pattern
- solve the problem: implements the time-advancing scheme and calls `getfem::standard_solve(...)`, which has already Netwon-like solver implemented and optimally chooses the best solver depending on the information given in the `assemble()` method.
- export the solution in `vtk` format.
- import and export the solution in `csv` format, mainly for inizialization and error computation purposes.


## Run the code
The user interfaces with the code just by modifying the getpot datafile in the `./examples/user/` directory. To execute the code the user just needs to run
```bash
./main <-v> <-m>
```
with `-v` if you want verbose output and `-m` if you want to create the mesh using Gmsh (in such case you need to have it installed in your system).
To check all the possible options, run
```bash
./main -h
```

Notice that if the mesh isn't generated using gmsh, the angle parameter will be ignored.
Depending on the method choosen, only the relative numerical parameters will be used.
As a reference manual for the code, check the "Implementation" section of the report.