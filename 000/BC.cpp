#include "BC.hpp"


namespace gf {

    BC::BC(Domain& domain) : M_domain(domain) {
        getfem::outer_faces_of_mesh(M_domain.getMesh(), border_faces);
    }

    BCNeu::BCNeu(const GetPot& datafile, bool isN)
    : M_BCType(BCType::Neumann)
    {
        // for (getfem::mr_visitor it(border_faces); it.finished(); ++it) {
        //     assert(it.is_face());
        //     base_node un = mesh.normal_of_face_of_convex(it.cv(), it.f());
        //     un /= gmm::vect_norm2(un);

        //     if (gmm::abs(un[0] - 1.0) < 1.0E-7)
        //     mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(), it.f());
        //     else if (gmm::abs(un[0] + 1.0) < 1.0E-7) 
        //     mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
            
        // }
        std::istringstream stream;
        if (isN) {
            std::string where = datafile("physics/regionNeuN", "10");
            M_expr = datafile("physics/bdLoadN", ...);
        }
        else {
            std::string where = datafile("physics/regionNeuT", "7");
            M_expr = datafile("physics/bdLoadT", ...);
        }
        
    }

    class BCbuilder {
        BC& M_BC;
        getfem::mesh_region
    }


}