#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include "Core.hpp"
#include "GetPot"

namespace gf {

    struct Domain {
        scalar_type Lx;
        scalar_type Ly;
        scalar_type Lz;
        size_type Nx;
        size_type Ny;
        size_type Nz;
    };
    struct Physics {
        scalar_type M_E;
        scalar_type M_nu;
        scalar_type M_lambda;
        scalar_type M_mu;
        base_small_vector M_gravity;
    };
    struct It {
        scalar_type maxIt;
        scalar_type rtol;
        scalar_type atol;
        /* other parameters ...*/
    };
    struct Nitsche {
        scalar_type theta;
        scalar_type gamma0;
    };
    struct Time {
        scalar_type t0;
        scalar_type tend;
        scalar_type dt;
    };

    struct Params {
        Domain domain;
        Physics phyisics;
        It it;
        Nitsche nitsche;
        Time time;

        Params(const GetPot&);
    };

} // namespace gf


#endif // _PARAMS_HPP_