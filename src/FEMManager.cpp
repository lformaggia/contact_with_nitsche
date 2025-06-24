#include "FEMManager.hpp"

namespace gf {

    FEMManager::FEMManager(const getfem::mesh& mesh, bool v)
        : M_mfU(mesh, mesh.dim()),
        M_mfStress(mesh),
        M_mfRhs(mesh),
        M_mfLM(mesh, mesh.dim()),
        verbose(v){
        }

    void FEMManager::setMeshFem(const Numerics &n){
        
        if (verbose) std::cout << "Setting Finite Element...";

        getfem::pfem pfU = getfem::fem_descriptor(n.FEMTypeDisplacement);
        getfem::pfem pfStress = getfem::fem_descriptor(n.FEMTypeStress);
        getfem::pfem pfRhs = getfem::fem_descriptor(n.FEMTypeRhs);
        getfem::pfem pfLM = getfem::fem_descriptor(n.FEMTypeLM); // Linear multipliers

        M_mfU.set_finite_element(pfU);
        M_mfStress.set_finite_element(pfStress);
        M_mfStress.set_qdim(3,3);
        M_mfRhs.set_finite_element(pfRhs);
        M_mfLM.set_finite_element(pfLM);

        if (verbose) std::cout << "done." << std::endl;

    }

} // namespace gf