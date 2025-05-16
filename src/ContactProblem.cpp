#include "ContactProblem.hpp"

namespace gf{

    ContactProblem::ContactProblem(const Mesh& m, const Params& p)
    : M_mesh(m),
    M_params(p),
    M_BC(m.get()),
    M_integrationMethod(m.get()),
    M_imData(M_integrationMethod)
    {

    }

    void ContactProblem::init() {

        M_BC.readBC(M_params.datafile, M_mesh.getBdRegions());

        M_FEM.setMeshFem(M_params.datafile, M_mesh.get(), M_mesh.getRegions());
                
        getfem::pintegration_method ppi = getfem::int_method_descriptor(M_params.numerics.integration);
        size_type N = M_params.domain.dim;
        M_integrationMethod.set_integration_method(M_mesh.get().convex_index(), ppi);
        M_imData.set_tensor_size(bgeot::multi_index(N,N));




    }


}