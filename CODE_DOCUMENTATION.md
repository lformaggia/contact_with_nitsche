# Contact with Nitsche - C++ Code Documentation

## Overview

This is a comprehensive C++ finite element framework for solving **contact mechanics problems with frictional interfaces** (fault modeling). The code implements multiple numerical methods for enforcing non-penetrability and Coulomb friction conditions on 3D domains with embedded faults. It was developed by S. Galati for the A.A. 2024-2025 academic year as part of a scientific computing project.

## Objectives

The primary objectives of this codebase are:

1. **Contact Mechanics Simulation**: Solve 3D contact problems on parallelepiped domains intersected by fault planes
2. **Multiple Enforcement Methods**: Implement and compare different approaches for contact and friction:
   - **Nitsche Method**: A consistent, non-conforming method for contact enforcement
   - **Penalty Method**: Approximate enforcement using penalty parameters
   - **Augmented Lagrangian**: Mixed formulation combining Lagrange multipliers with penalty terms
3. **Fault Modeling**: Handle slip and contact conditions along geological fault interfaces
4. **Mesh Generation**: Support both built-in structured meshes and external Gmsh integration for complex geometries
5. **Scientific Computing**: Demonstrate advanced finite element techniques using the GetFEM++ library

## Architecture and Design

The codebase follows a **modular, object-oriented design** with clear separation of concerns and extensive use of design patterns:

### Core Design Patterns

- **Strategy Pattern**: `ContactEnforcementStrategy` hierarchy allows runtime selection of contact methods
- **Builder Pattern**: `MeshBuilderStrategy` enables different mesh creation approaches
- **Factory Pattern**: Dynamic creation of contact enforcement and mesh building strategies
- **RAII Pattern**: Extensive use of smart pointers and automatic resource management

### Key Components

#### 1. Parameter Management (`Params.hpp/cpp`)

```cpp
struct Params {
    Domain domain;      // Geometry and discretization
    Physics physics;    // Material properties and loading
    Contact contact;    // Contact method and parameters
    Time time;         // Temporal discretization
    Numerics numerics; // FE spaces and integration
    // ... flags and options
};
```

**Purpose**: Centralized configuration management using GetPot library for parameter parsing.

**Key Features**:
- Reads from `data.pot` configuration files
- Command-line argument processing
- Automatic parameter validation and derived quantities calculation
- Support for different mesh types (tetrahedra vs hexahedra)

#### 2. Mesh Management (`Mesh.hpp/cpp`, `MeshBuilder.hpp/cpp`)

```cpp
class Mesh {
    std::unique_ptr<MeshBuilderStrategy> M_meshBuilder;
    getfem::mesh M_mesh;
    // ...
};

class MeshBuilderStrategy {
    virtual void buildMesh(getfem::mesh&) const = 0;
    virtual void initRegions(getfem::mesh&) const = 0;
};
```

**Purpose**: Flexible mesh creation and region management.

**Strategies**:
- `BuiltInBuilder`: Creates structured meshes using GetFEM++ primitives
- `GmshBuilder`: Imports meshes from Gmsh `.msh` files, supports complex geometries

**Region Management**:
- `BulkLeft`/`BulkRight`: Separate domains on each side of the fault
- `Fault`: Interface region where contact conditions are enforced
- Boundary regions for applying boundary conditions

#### 3. Finite Element Management (`FEMManager.hpp/cpp`)

```cpp
class FEMManager {
    getfem::mesh_fem M_mfU;      // Displacement field
    getfem::mesh_fem M_mfStress; // Stress field  
    getfem::mesh_fem M_mfRhs;    // Right-hand side
    getfem::mesh_fem M_mfLM;     // Lagrange multipliers
};
```

**Purpose**: Manages multiple finite element spaces for different field variables.

**Features**:
- Supports different polynomial orders for each field
- Handles both P1/Q1 (linear) and higher-order elements
- Automatic setup based on mesh type (tetrahedra/hexahedra)

#### 4. Boundary Condition Framework (`BCHandler.hpp/cpp`, `BC.hpp`)

```cpp
class BC {
    getfem::mesh_region M_region;
    VectorFunctionType M_function;  // Mathematical expression
    BCType M_BCtype;               // Dirichlet/Neumann/Mixed
};

class BCHandler {
    std::unordered_map<BCType, std::vector<std::unique_ptr<BC>>> M_BCList;
    muParserXInterface M_parser;   // Expression parsing
};
```

**Purpose**: Flexible boundary condition specification and management.

**Supported Types**:
- **Dirichlet**: Fixed displacement conditions
- **Neumann**: Applied traction/force conditions  
- **Mixed**: Normal displacement constraint (allows tangential slip)

**Key Features**:
- Mathematical expression parsing using muParserX
- Time-dependent boundary conditions
- Automatic region-based application

#### 5. Contact Enforcement Strategies (`ContactEnforcementStrategy.hpp/cpp`)

The heart of the contact mechanics implementation:

```cpp
class ContactEnforcementStrategy {
    virtual void enforce(getfem::model& md, const getfem::mesh_im& im, bool verbose) const = 0;
};
```

##### Nitsche Method (`NitscheContactEnforcement`)
- **Mathematical Foundation**: Consistent weak formulation without Lagrange multipliers
- **Parameters**: 
  - `θ` (theta): Symmetry parameter (-1: skew-symmetric, 0: nonsymmetric, 1: symmetric)
  - `γN` (gammaN): Penalty-like parameter for stability
- **Advantages**: Optimal convergence rates, no additional unknowns
- **Implementation**: Uses GetFEM++'s Generic Weak Form Language (GWFL)

##### Penalty Method (`PenaltyContactEnforcement`)
- **Mathematical Foundation**: Approximate enforcement using large penalty parameters
- **Parameters**: `γP` (gammaP): Penalty parameter
- **Advantages**: Simple implementation, no Lagrange multipliers
- **Disadvantages**: Conditioning issues, approximate satisfaction of constraints

##### Augmented Lagrangian (`AugmentedLagrangianContactEnforcement`)
- **Mathematical Foundation**: Mixed formulation combining Lagrange multipliers with penalty
- **Parameters**: `γL` (gammaL): Augmentation parameter  
- **Advantages**: Exact constraint enforcement, better conditioning than pure penalty
- **Implementation**: Requires additional mesh finite element space for multipliers

#### 6. Main Problem Class (`ContactProblem.hpp/cpp`)

```cpp
class ContactProblem {
    const Mesh& M_mesh;
    const Params& M_params;
    BCHandler M_BC;
    FEMManager M_FEM;
    std::unique_ptr<ContactEnforcementStrategy> M_contactEnforcement;
    getfem::model M_model;  // GetFEM++ model
};
```

**Purpose**: Orchestrates the entire simulation workflow.

**Workflow**:
1. **`init()`**: Setup FE spaces, integration methods, boundary conditions
2. **`assemble()`**: Define the weak formulation using GWFL
3. **`solve()`**: Time-stepping loop with Newton solver

## Mathematical Formulation

### Elasticity Problem

The code solves the **linear elasticity equations** on domains $\Omega_L$ and $\Omega_R$ separated by fault $\Gamma$:

$$
\begin{align}
-\nabla \cdot \sigma(\mathbf{u}) &= \mathbf{f} \quad \text{in } \Omega_L \cup \Omega_R \\
\sigma(\mathbf{u}) &= \lambda \text{tr}(\varepsilon(\mathbf{u})) \mathbf{I} + 2\mu \varepsilon(\mathbf{u}) \\
\varepsilon(\mathbf{u}) &= \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T)
\end{align}
$$

### Contact and Friction Conditions

On the fault interface $\Gamma$:

1. **Non-penetration**: $[u_n] \leq 0$ (gap condition)
2. **Complementarity**: $\sigma_n \leq 0$, $\sigma_n [u_n] = 0$
3. **Coulomb friction**: $|\sigma_t| \leq \mu_f |\sigma_n|$

Where $[u] = u_L - u_R$ denotes the jump across the interface.

### Nitsche Formulation

The Nitsche method adds the following terms to the weak formulation:

$$
\int_\Gamma \left( \{[\sigma(\mathbf{u})]\} \cdot [\mathbf{v}] - \theta [\mathbf{u}] \cdot \{[\sigma(\mathbf{v})]\} + \frac{\gamma}{h} [\mathbf{u}] \cdot [\mathbf{v}] \right) d\Gamma
$$

Plus contact and friction nonlinearities.

## Implementation Details

### GetFEM++ Integration

The code extensively uses **GetFEM++** features:

- **Generic Weak Form Language (GWFL)**: Symbolic definition of integrals
- **Model Framework**: Automatic linearization and assembly
- **Built-in Solvers**: Newton methods with line search
- **Integration Methods**: High-order quadrature rules

### Key Macros and Expressions

The implementation defines numerous GWFL macros for contact terms:

```cpp
// Jump operators
M_model.add_macro("u_jump", "(uL - Interpolate(uR,neighbor_element))");
M_model.add_macro("un_jump", "u_jump . n");

// Stress computations  
M_model.add_macro("stressL", "(lambda*Trace(Grad_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_uL))");

// Contact conditions
M_model.add_macro("Pn_u", "(gammaN * un_jump - sig_u_nL)");
```

### Friction Modeling

Coulomb friction is implemented using **projection operators**:

```cpp
M_model.add_macro("Sh", "mu_fric * pos_part(Pn_u)");  // Friction threshold
M_model.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
```

## Technical Features

### Mesh Generation
- **Built-in**: Structured hexahedral/tetrahedral meshes
- **Gmsh Integration**: Supports complex geometries, inclined faults
- **Region Management**: Automatic fault detection and boundary labeling

### Time Integration
- **Quasi-static**: Load stepping for nonlinear contact problems
- **Adaptive Loading**: Time-dependent boundary conditions and body forces

### Solution Export
- **VTK Format**: Paraview-compatible visualization with stress fields
- **CSV Export**: Numerical data for post-processing and validation
- **Displacement and Stress**: Full field output with proper fault visualization

### Error Estimation
- **Convergence Testing**: Built-in mesh refinement and error computation
- **Reference Solutions**: CSV import/export for validation

## Dependencies and Build System

### Required Libraries
- **GetFEM++**: Core finite element library
- **GMM++**: Linear algebra (part of GetFEM++)
- **MUMPS**: Direct sparse solver
- **SuperLU**: Alternative sparse solver  
- **LAPACK/BLAS**: Dense linear algebra
- **Qhull**: Computational geometry
- **muParserX**: Mathematical expression parsing

### CMake Build System
```cmake
# Static library creation
add_library(mycontactlib STATIC ${CORE_SOURCES})

# External dependencies
ExternalProject_Add(muparserx_ext ...)

# Example automation
foreach(EXAMPLE ${EXAMPLE_DIRS})
    add_test(NAME run_${EXAMPLE} ...)
endforeach()
```

**Features**:
- External project management for muParserX
- Automatic test generation for examples
- Cross-platform compatibility (Unix symlinks/Windows copy)
- Doxygen documentation generation

## Usage and Examples

### Command Line Interface
```bash
./main [-v] [-m] [-f datafile] [-i] [-r] [-t]
```

**Options**:
- `-v`: Verbose output
- `-m`: Use Gmsh for mesh generation  
- `-f`: Specify data file (default: `data.pot`)
- `-i`: Initialize from previous solution
- `-r`: Generate refined reference solution
- `-t`: Test mode (compute errors)

### Configuration File Format (`data.pot`)

```bash
# Domain parameters
domain/Lx = 2.0          # Domain length in x
domain/Ly = 1.0          # Domain length in y  
domain/Lz = 1.0          # Domain length in z
domain/h = 0.1           # Characteristic mesh size
domain/meshType = "GT_PK(3,1)"  # Tetrahedra or "GT_QK(3,1)" for hexahedra

# Physics parameters
physics/E = 1e6          # Young's modulus
physics/nu = 0.3         # Poisson's ratio
physics/mu_friction = 0.5 # Friction coefficient
physics/bulkLoad = "[0., 0., -9.81]"  # Body force (gravity)

# Contact method
contact/method = "nitsche"  # "penalty" or "augLM"  
contact/theta = 0.0         # Nitsche parameter
contact/gammaN = 10.0       # Penalty parameter (scaled by E/h)

# Time parameters
time/t0 = 0.0
time/tend = 1.0
time/dt = 0.1

# Boundary conditions
regionDirichlet = [1 2]     # Boundary region IDs
bdDirichlet1 = "{0,0,0}"    # Fixed displacement
regionNeumann = [5 6]       # Traction regions  
bdNeumann1 = "{0,0,-1000}"  # Applied traction
```

## Scientific Applications

This code framework is particularly well-suited for:

### Geomechanics
- **Fault modeling**: Earthquake mechanics, fault slip analysis
- **Rock mechanics**: Fracture propagation, joint behavior
- **Reservoir geomechanics**: Hydraulic fracturing, subsidence

### Engineering Applications  
- **Contact mechanics**: Interface problems in mechanical engineering
- **Fracture mechanics**: Crack propagation with contact
- **Multiphysics**: Foundation for coupled problems (hydro-mechanical, thermo-mechanical)

### Numerical Methods Research
- **Contact algorithms**: Comparison of enforcement methods
- **Finite element methods**: Non-conforming methods, mixed formulations
- **Computational mechanics**: Large-scale 3D simulations

## Code Quality and Documentation

### Documentation Standards
- **Doxygen**: Comprehensive API documentation
- **Inline Comments**: Mathematical formulations and algorithm explanations
- **README**: Build instructions and usage examples

### Software Engineering
- **RAII**: Automatic memory management with smart pointers
- **Exception Safety**: Proper error handling throughout
- **Const Correctness**: Immutable interfaces where appropriate
- **Modular Design**: Clear separation of concerns

### Testing and Validation
- **CTest Integration**: Automated example testing
- **Convergence Studies**: Built-in error computation capabilities
- **Multiple Examples**: Various contact scenarios for validation

This codebase represents a sophisticated implementation of contact mechanics algorithms, demonstrating advanced C++ programming techniques, numerical methods, and scientific computing best practices. It serves both as a research tool and an educational example of finite element software development.