#ifndef _FEM_MANAGER_HPP_
#define _FEM_MANAGER_HPP_

#include "Core.hpp"
#include "Params.hpp"

namespace gf {

    /**
     * @brief Class that manages the finite element method (FEM) spaces for the problem.
     * It creates and stores mesh_fem objects for displacement, stress, right-hand side, and Lagrange multipliers.
     */
    class FEMManager {
    public:
        FEMManager(const getfem::mesh&);
    
        /**
         * @brief Set the finite element method (FEM) spaces based on the provided numerics and mesh.
         * @param n The Numerics object containing FEM type information.
         * @param m The mesh to set the FEM spaces for.
         */
        void setMeshFem(const Numerics& n, const getfem::mesh&);

        getfem::mesh_fem& mf_u() { return M_mfU; }
        
        getfem::mesh_fem& mf_stress() { return M_mfStress; }

        getfem::mesh_fem& mf_rhs() { return M_mfRhs; }

        getfem::mesh_fem& mf_LMn() { return M_mfLMn; }
        
        getfem::mesh_fem& mf_LMt() { return M_mfLMt; }

    private:
        getfem::mesh_fem M_mfU;
        getfem::mesh_fem M_mfStress;
        getfem::mesh_fem M_mfRhs;
        getfem::mesh_fem M_mfLMn;
        getfem::mesh_fem M_mfLMt;

    };

}



#endif // _FEM_MANAGER_HPP_