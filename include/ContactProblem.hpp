#ifndef _CONTACT_PROBLEM_HPP_
#define _CONTACT_PROBLEM_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "Mesh.hpp"
#include "BCHandler.hpp"
#include "FEMManager.hpp"
#include "ContactEnforcementStrategy.hpp"
#include "GetPot"

namespace gf{

    /** \todo: update
     * @brief The main class that orchestrates the problem for contact mechanics.
     * It stores the GetPot and the datafile, which needs to be shared between the different classes.
     * It allows to perfom 3d simulation on a 3D parallelepiped, with a structured conforming mesh, cut by
     * a plane at x=0 with normal aligned to the x-axis.
     */
    class ContactProblem {
    public:
    
        /**
         * @brief Constructor taking the data filename
         * 
         */
        ContactProblem (const Mesh&, const Params&);
        
        /**
         * @brief Initialize the problem
         * 1. import BCs
         * 2. set FE spaces
         * 3. set integration method
         */
        void init();

        /**
         * @brief Assemble the linear system based on the GWFL (Generic Weak Form Language)
         */
        void assemble();

        /**
         * @brief Solve the linear system
         */
        void solve();

        
    private:

        const Params& M_params;
        const Mesh& M_mesh; ///< The mesh
        BCHandler M_BC; ///< Class that stores BC information
        FEMManager M_FEM; ///< Stores the getfem::mesh_fem objects
        std::unique_ptr<ContactEnforcementStrategy> M_contactEnforcement; ///< Strategy for contact enforcement
        getfem::mesh_im M_integrationMethod; ///< Integration methods
        getfem::im_data M_imData;
        getfem::model M_model;

        /**
         * @brief Export vtk results for visualization
         */
        void exportVtk(size_type i);



    };


}

#endif // _CONTACT_PROBLEM_HPP_