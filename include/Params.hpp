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
        std::string meshType;
    };
    struct Physics {
        scalar_type M_E0;
        scalar_type M_nu;
        scalar_type M_lambda;
        scalar_type M_mu;
        base_small_vector M_gravity;
    };
    struct It {
        scalar_type maxIt;
        scalar_type rtol;
        scalar_type atol;
        /** \todo:
         * add other Newton Params, if needed
         */
    };
    struct Nitsche {
        scalar_type theta;
        scalar_type gamma0;
        /** \todo:
         * add tolerances for Nitsche
         */
    };
    struct Numerics {
        std::string integration;
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
        Nitsche nitsche;
        Time time;
        std::string meshFile;
        Numerics numerics;
        bool verbose;
        
        Params(const std::string&, const std::string&, bool);
        
    };

    
    std::ostream& operator<<(std::ostream&, const Params&);

} // namespace gf


#endif // _PARAMS_HPP_