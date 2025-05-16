#include "FEMManager.hpp"

bool DEBUGFEM = true;
namespace gf {

    void FEMManager::setMeshFem(const GetPot& datafile, const getfem::mesh& mesh, const RegionMapType& regions){
        std::string FEMTypeDisp = datafile("numerics/FEMTypeDisplacement", "FEM_QK(3,1)");
        std::string FEMTypeStress = datafile("numerics/FEMTypeStress", "FEM_QK(3,1)");
        std::string FEMTypeRhs = datafile("numerics/FEMTypeRhs", "FEM_QK(3,1)");
        std::string FEMTypeCoeff = datafile("numerics/FEMTypeCoeff", "FEM_QK(3,0)");

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
        if (DEBUGFEM){
            std::clog << "Setting Finite Element... ";
        }
        M_mfU1.set_finite_element(regions.at("BulkLeft")->index(), pfU);
        M_mfU2.set_finite_element(regions.at("BulkRight")->index(), pfU);
        M_mfStress1.set_finite_element(regions.at("BulkLeft")->index(), pfStress);
        M_mfStress2.set_finite_element(regions.at("BulkRight")->index(), pfStress);
        M_mfRhs.set_finite_element(pfRhs);
        M_mfCoeff.set_finite_element(pfCoeff);

        if (DEBUGFEM){
            std::clog << "done." << std::endl;
        }

    }

} // namespace gf