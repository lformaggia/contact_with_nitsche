#ifndef _DOMAIN_HPP_
#define _DOMAIN_HPP_

#include "Core.hpp"
#include "DomainView.hpp"
#include <unordered_map>


namespace gf
{

    class MeshBuilderStrategy;
    class GmshBuilder;
    class BuiltInBuilder;

    /**
     * @brief Class that stores all the domain information
     * Key information:
     * - it owns the mesh object and delegates its building to a MeshBuilderStrategy class
     * - supports both 2D/3D meshes
     * - stores the physical parameters (lamé coefficients)
     * - 
     */
    class Domain {
    private:
        friend class BuiltInBuilder;
        friend class GmshBuilder;
        
    public:
        /**
         * @brief Constructor to createe the domain from the datafile
         * @param datafile the GetPot object containing all the input data and parameters
         * Depending on the value of meshFile, it delegates the construction to meshBuilder
         */
        Domain(const GetPot& datafile);

        /**
         * @brief Returns a const reference to the mesh object (read-only)
         */
        const getfem::mesh& getMesh() const { return M_mesh; }

        /**
         * @brief Getter for mu
         */
        inline scalar_type mu() const { return M_mu; }
        
        /**
         * @brief Getter for lambda
         */
        inline scalar_type lambda() const { return M_lambda; }

        
        /**
         * @brief Export the mesh for visualization
         * @param Filename The output filename
         */
        void exportMesh(const std::string& filename) const;
        
    private:
        std::unordered_map<std::string, std::unique_ptr<DomainView>> M_regions; ///< Views on the different regions
        std::unique_ptr<MeshBuilderStrategy> M_meshBuilder; ///< MeshBuilder
        getfem::mesh M_mesh; ///< The mesh object
        std::string M_meshType; ///< mesh type (default: linear tetrahedra)
        scalar_type M_mu; ///< First lamé parameter
        scalar_type M_lambda; ///< Second lamé parameter
    };

}






#endif // _DOMAIN_HPP_