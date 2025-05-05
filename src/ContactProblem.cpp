#include "ContactProblem.hpp"

namespace gf{

    ContactProblem::ContactProblem (const std::string& filename)
    : M_datafile(filename.data()){}

    void ContactProblem::init() {
        if (true /*MODIFY THIS*/)
            M_meshBuilder = std::make_unique<BuiltInBuilder>(M_datafile);
        M_meshBuilder->buildMesh(M_mesh);
        M_meshBuilder->initRegions(M_mesh, M_regions);

        BCHandler.setBC(M_datafile);
        FEManager.setFEM();
    }

    

}