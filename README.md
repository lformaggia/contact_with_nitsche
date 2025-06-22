# A NITSCHE-BASED METHOD FOR SOLVING A PROBLEM IN CONTACT MECHANICS - Sequential version
## S. Galati
### A.A. 2024-2025

The project aims to solve a contact mechanics problem using a Nitsche-based method for imposing both the non-penetrability and the Coulomb friction condition at a fault cutting a 3-dimensional bulk domain.
The FEM framework is given by GetFEM++, a very powerful library for solving many FEM problems. In order to compile and execute the code, the user will need the following libraries installed:
- getfem
- mumps (sequential version): for the linear algebra procedures

The mesh can be built either internally to getfem or by using gmsh and giving the option `-m` at runtime (in such a case, gmsh needs to be installed). The user has just to modify the data.pot file properly.

The `src` and `include` directories contain the code, `external` just contains a copy of muparserx (which is used to read functions form the .pot file), while the `examples` directory contains some tests and the folder where the user should work.

### Source code
The main implemented classes are:
`Params`:
A struct holding all the parameters contained in `data.pot`.
- It implements a function for reading the parameters from a GetPot datafile, including some options passed at runtime.
- Overloads the streaming operator.

`Mesh` and `MeshBuilderStrategy`: the first is just a wrapper for the `getfem::mesh` object, which is built using the Builder in `MeshBuilder.hpp`. It allows direct access to the underlying regions definied within the mesh (Bulks, Fault and Boundary regions).

`FEMManager`: to manage the `getfem::mesh_fem` objects.

`BCHandler` and `BC`: to handle the boundary conditions. The supported BC types are Dirichlet, Neumann and normal Dirichlet. In the data.pot:
- region< type > = [1 4 5]: Dirichlet is imposed on boundary regions with tags 1,4 and 5
- bd< type >< i > = {x[1]+sqrt(x[0]),x[2],3.5*t^2}: the < type > to be imposed on the i-th region given in region< type >

`ContactProblem` and `ContactEnforcementStrategy`: the main class handling the problem. It implements methods to:
- initialize the problem: set BCs, set FEM, set integration method
- prepare the assemble: use GFWL of getfem to instruct the solver on how to assemble the matrices at each iteration, choosing the method through a strategy pattern
- solve the problem: implements the time-advancing scheme and calls `getfem::standard_solve(...)`, which has already Netwon-like solver implemented and optimally chooses the best solver depending on the information given in the `assemble()` method.
- export the solution in `vtk` format.


## Notes for compiling the code
Before compiling the code, make sure you have the getfem (sequential version, static library) and mumps (sequential version, static library) installed. In order to avoid conflicts between uncompatible versions of some libraries, ensure that you have no modules loaded, otherwise just run
```bash
module purge
```

TODO: Give information on how to install those libraries


To compile the code, in the current repository run
```bash
cmake -B build
```
with possible options (all defaulted to OFF):
- `-DEXAMPLES_VERBOSE=ON`: to be verbose when testing the examples
- `-DEXAMPLES_USE_GMSH=ON`: to run the examples using a mesh generated via gmsh

and then run:
```bash
cmake --build build -j<nprocs>
```

This will build muparserx library and install a shared version that the code links againts, build a static library `libmycontactlib.a`, which is used by the executable, and create the executable. Notice that the executable is stored in the `./build/bin/` directory and then a symbolic link (or a copy if the system is on non-Unix type) is created in each `./examples/<example>` folder.

If you want to run the examples for testing the code, optionally run
```bash
cd build && ctest --output-on-failure
```
This will run all the predefined examples, printing the ouput of execution in a `execution_log.txt` file and the vtk files for visualization in its `output/` subfolder

To clean the code, the following target has been defined (from the root of the project):
```bash
cmake --build build --target clean-all
```
This will remove the `build/` and `lib/` directory, together with the output of the execution of the examples if `ctest` has been runned.


## Run the code
The user interfaces with the code just by modifying the getpot datafile in the `./examples/user/` directory. To execute the code the user just needs to run
```bash
./main
```
with possible options:
  - `-v`: be verbose
  - `-m`: generate the mesh with gmsh

Notice that if the mesh isn't generated using gmsh, the angle parameter will be ignored.
Depending on the method choosen, only the relative numerical parameters will be used.