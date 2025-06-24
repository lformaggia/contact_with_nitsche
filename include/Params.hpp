#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include "Core.hpp"
#include "GetPot"

namespace gf {

    /**
     * @brief A structure to hold the domain parameters of the simulation
     * It contains the dimensions, lengths, number of elements in each direction,
     * mesh type, and other related parameters.
     */
    struct Domain {
        size_type dim;
        scalar_type Lx;
        scalar_type Ly;
        scalar_type Lz;
        size_type Nx;
        size_type Ny;
        size_type Nz;
        scalar_type h;
        scalar_type angle;
        std::string meshType;
    };

    /**
     * @brief A structure to hold the physics parameters of the simulation
     * It contains the material properties such as Young's modulus, Poisson's ratio,
     * and the gravitational acceleration.
     */
    struct Physics {
        scalar_type M_E0;
        scalar_type M_nu;
        scalar_type M_lambda;
        scalar_type M_mu;
        base_small_vector M_gravity;
        scalar_type M_mu_friction;
    };

    /**
     * @brief A structure to hold the parameters related to the Newton method
     * It contains the maximum number of iterations and the tolerance for convergence.
     * Additional parameters can be added as needed.
     */
    struct It {
        scalar_type maxIt;
        scalar_type tol;
        /** Possible future extension:
         * add other Newton Params, if needed
         */
    };

    /**
     * @brief A structure to hold the contact parameters of the simulation
     * It contains the method used for contact, and the parameters related to that method.
     * The method can be "penalty", "nitsche", or "augLM".
     */
    struct Contact {
        std::string method;
        scalar_type theta;
        scalar_type gammaN;
        scalar_type gammaP;
        scalar_type gammaL;
    };

    /**
     * @brief A structure to hold the numerical parameters of the simulation
     * It contains the integration method, the finite element types for displacement, stress, rhs, and LM.
     */
    struct Numerics {
        std::string integration;
        std::string FEMTypeDisplacement;
        std::string FEMTypeStress;
        std::string FEMTypeRhs;
        std::string FEMTypeLM;
    };

    /**
     * @brief A structure to hold the time parameters of the simulation
     * It contains the initial time, final time, and time step.
     * It can be extended to include backtracking parameters if needed.
     */
    struct Time {
        scalar_type t0;
        scalar_type tend;
        scalar_type dt;
        /** Possible future extension:
         * add backtracking parameters, if needed
         */
    };

    /**
     * @brief A structure to hold all the parameters of the simulation
     * It is initialized from the command line arguments and the data file
     * It contains all the parameters needed for the simulation, including
     * domain, physics, time, numerics, etc.
     */
    struct Params {
        GetPot datafile;
        Domain domain;
        Physics physics;
        It it;
        Contact contact;
        Time time;
        Numerics numerics;
        bool verbose = false;
        bool gmsh = false;
        
        Params(int argc, char* argv[]);
        
    };

    /**
     * @brief Overloaded output operator for Params
     */
    std::ostream& operator<<(std::ostream&, const Params&);

} // namespace gf


#endif // _PARAMS_HPP_