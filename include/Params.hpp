#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include "Core.hpp"
#include "GetPot"

namespace gf {

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
    struct Physics {
        scalar_type M_E0;
        scalar_type M_nu;
        scalar_type M_lambda;
        scalar_type M_mu;
        base_small_vector M_gravity;
        scalar_type M_mu_friction;
    };
    struct It {
        scalar_type maxIt;
        scalar_type rtol;
        scalar_type atol;
        /** \todo:
         * add other Newton Params, if needed
         */
    };
    struct Contact {
        std::string method;
        scalar_type theta;
        scalar_type gammaN;
        scalar_type gammaP;
        scalar_type gammaL;
    };
    struct Numerics {
        std::string integration;
        std::string FEMTypeDisplacement;
        std::string FEMTypeStress;
        std::string FEMTypeRhs;
        std::string FEMTypeLM;
    }; 
    struct Time {
        scalar_type t0;
        scalar_type tend;
        scalar_type dt;
        /** \todo:
         * add backtracking parameters, if needed
         * */
    };

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

    
    std::ostream& operator<<(std::ostream&, const Params&);

} // namespace gf


#endif // _PARAMS_HPP_