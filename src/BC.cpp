#include "BC.hpp"

namespace gf {

    BC::BC(const getfem::mesh_region& region, size_type regionID, VectorFunctionType&& func)
    : M_region(region), M_ID(regionID), M_function(std::move(func))
    {

    }



    BCNeu::BCNeu(const GetPot& datafile, bool isN)
    : M_BCType(BCType::Neumann)
    {
        if (isN) {
            std::string where = datafile("physics/regionNeuN", "10");
            M_expr = datafile("physics/bdLoadN", ...);
        }
        else {
            std::string where = datafile("physics/regionNeuT", "7");
            M_expr = datafile("physics/bdLoadT", ...);
        }
        
    }





} // namesace gf