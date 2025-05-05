#include "ContactProblem.hpp"

namespace gf{

    ContactProblem::ContactProblem (const std::string& filename)
    : M_datafile(filename.c_str()), M_params(M_datafile){}

    void ContactProblem::init() {
        if (true /*MODIFY THIS*/)
            M_meshBuilder = std::make_unique<BuiltInBuilder>(M_datafile);
        M_meshBuilder->buildMesh(M_mesh);
        M_meshBuilder->initRegions(M_mesh, M_regions);

        BCHandler.readBC(M_datafile, M_params.domain);
        FEManager.setFEM();
        
    }

    

}