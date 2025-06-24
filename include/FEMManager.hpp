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

        /**
         * @brief Constructor that initializes the FEMManager with a mesh.
         * @param m The mesh to associate with this FEMManager.
         * @param v Verbosity flag to control logging output.
         */
        FEMManager(const getfem::mesh& m, bool v = false);
    
        /**
         * @brief Set the finite element method (FEM) spaces based on the provided numerics and mesh.
         * @param n The Numerics object containing FEM type information.
         */
        void setMeshFem(const Numerics& n);

        /**
         * \defgroup MeshFemGetters Getters for mesh_fem
         * \brief Getters that return references to various mesh_fem objects.
         *
         * These functions provide access to internal mesh_fem representations
         * used for displacement, stress, RHS, and Lagrange multipliers.
         */
        ///@{
        getfem::mesh_fem& mf_u() { return M_mfU; }
        
        getfem::mesh_fem& mf_stress() { return M_mfStress; }

        getfem::mesh_fem& mf_rhs() { return M_mfRhs; }

        getfem::mesh_fem& mf_LM() { return M_mfLM; }
        ///@}

    private:
        getfem::mesh_fem M_mfU; ///< Mesh finite element method for displacement
        getfem::mesh_fem M_mfStress; ///< Mesh finite element method for stress
        getfem::mesh_fem M_mfRhs; ///< Mesh finite element method for right-hand side
        getfem::mesh_fem M_mfLM; ///< Mesh finite element method for Lagrange multipliers
        bool verbose; ///< Verbosity flag for logging

    };

}



#endif // _FEM_MANAGER_HPP_