#include "BC.hpp"

namespace gf {

    BC::BC(const Boundary& bdRegion, size_type regionID, VectorFunctionType func, BCType bctype)
    : M_region(bdRegion), M_ID(regionID), M_function(func), M_BCtype(bctype)
    {
    }

    // BCNeu::BCNeu(const GetPot& datafile, bool isN)
    // : 
    // {
    //     if (isN) {
    //         std::string where = datafile("physics/regionNeuN", "10");
    //         M_expr = datafile("physics/bdLoadN", ...);
    //     }
    //     else {
    //         std::string where = datafile("physics/regionNeuT", "7");
    //         M_expr = datafile("physics/bdLoadT", ...);
    //     }
        
    // }


} // namesace gf