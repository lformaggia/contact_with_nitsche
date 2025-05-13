#include "FEMManager.hpp"

namespace gf {

    void FEMManager::setMeshFem(const GetPot& datafile, const getfem::mesh& mesh, const RegionMapType& regions){
        std::string FEMTypeDisp = datafile("FEMTypeDisplacement", "FEM_PK(3,1)");
        std::string FEMTypeStress = datafile("FEMTypeStress", "FEM_PK(3,1)");
        std::string FEMTypeRhs = datafile("FEMTypeRhs", "FEM_PK(3,1)");
        std::string FEMTypeCoeff = datafile("FEMTypeCoeff", "FEM_PK(3,0)");

        getfem::pfem pfU = getfem::fem_descriptor(FEMTypeDisp);
        getfem::pfem pfStress = getfem::fem_descriptor(FEMTypeStress);
        getfem::pfem pfRhs = getfem::fem_descriptor(FEMTypeRhs);
        getfem::pfem pfCoeff = getfem::fem_descriptor(FEMTypeCoeff);

        M_mfU1.init_with_mesh(mesh, 3);
        M_mfU2.init_with_mesh(mesh, 3);
        M_mfStress1.init_with_mesh(mesh,1);
        M_mfStress1.set_qdim(3,3);
        M_mfStress2.init_with_mesh(mesh, 1);
        M_mfStress2.set_qdim(3,3);
        M_mfRhs.init_with_mesh(mesh, 3);
        M_mfCoeff.init_with_mesh(mesh, 1);
        
        M_mfU1.set_finite_element(regions.at("Bulk1")->index(), pfU);
        M_mfU2.set_finite_element(regions.at("Bulk2")->index(), pfU);
        M_mfStress1.set_finite_element(regions.at("Bulk1")->index(), pfStress);
        M_mfStress2.set_finite_element(regions.at("Bulk2")->index(), pfStress);
        M_mfRhs.set_finite_element(pfRhs);
        M_mfCoeff.set_finite_element(pfCoeff);

    }

} // namespace gf