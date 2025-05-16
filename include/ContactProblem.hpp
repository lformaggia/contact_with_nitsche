#ifndef _CONTACT_PROBLEM_HPP_
#define _CONTACT_PROBLEM_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "MeshBuilder.hpp"
#include "MeshRegion.hpp"
#include "BCHandler.hpp"
#include "FEMManager.hpp"
#include "GetPot"

namespace gf{

    /**
     * @brief The main class that orchestrates the problem for contact mechanics.
     * It stores the GetPot and the datafile, which needs to be shared between the different classes.
     * It allows to perfom 3d simulation on a 3D parallelepiped, with a structured conforming mesh, cut by
     * a plane at x=0 with normal aligned to the x-axis.
     */
    class ContactProblem {
    public:
    
        /**
         * @brief Constructor taking the data filename
         * Initializes the stored GetPot object, which is needed by most of the methods
         */
        ContactProblem (const std::string& filename, const std::string& meshfile, bool verbose);
        
        /**
         * @brief Initialize the problem
         * 1. build the mesh
         * 2. build the mesh_regions
         * 3. import BCs
         * 4. set FE spaces
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

        /**
         * @brief Export vtk results for visualization
         */
        void exportResults() const;
        
    private:

        GetPot M_datafile;
        Params M_params;
        std::unique_ptr<MeshBuilderStrategy> M_meshBuilder; ///< MeshBuilder
        getfem::mesh M_mesh; ///< The mesh
        RegionMapType M_regions; ///< The mesh regions (BulkLeft, BulkRight, Fault)
        std::unique_ptr<BCHandler> M_BC; ///< Class that stores BC information
        FEMManager M_FEM; ///< Stores the getfem::mesh_fem objects
        getfem::mesh_im  M_IntegrationMethod; ///< Integration methods
        getfem::im_data M_imData;
        // TimeManager M_time; ///< Struct that keeps information for evolutionary problems
        // SolType M_U; ///< Vector to store the solution 
        // StressType M_stress; ///< Vector that stores the stress information

        /*...*/

    };

}

#endif // _CONTACT_PROBLEM_HPP_