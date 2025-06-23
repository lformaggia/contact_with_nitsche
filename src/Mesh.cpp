#include "Mesh.hpp"

namespace gf {

    Mesh::Mesh(const Params& p)
    : M_params(p)
    {
        if (M_params.gmsh)
            M_meshBuilder = std::make_unique<GmshBuilder>(M_params.domain);
        else
            M_meshBuilder = std::make_unique<BuiltInBuilder>(M_params.domain);
        
        M_meshBuilder->construct(M_mesh);

        if (M_params.verbose)
            M_mesh.write_to_file("mesh_getfem_info.log");

    }


} // namespace gf