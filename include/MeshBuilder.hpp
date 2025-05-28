
#ifndef _MESH_BUILDER_HPP_
#define _MESH_BUILDER_HPP_

#include "Core.hpp"
#include "Params.hpp"

namespace gf {

    class MeshRegion;

    /**
     * Builder class for building the mesh
     */
    class MeshBuilderStrategy{
    public:
        /**
         * @brief Null method to build the mesh, overridden in the hierarchy
         * They get a reference to the mesh object, modifying it
         */
        virtual void buildMesh(getfem::mesh&) const = 0;

        virtual void initRegions(getfem::mesh&) const = 0;

        ~MeshBuilderStrategy() = default;
    };


    class BuiltInBuilder : public MeshBuilderStrategy {
    private:
        Domain M_domain; ///< Domain parameters

    public:
        BuiltInBuilder(const Domain& d):
        M_domain(d){}

        /**
         * @brief Builds the mesh with the domain information
         */
        void buildMesh(getfem::mesh&) const override;

        /**
         * @brief Initializes regions for the domain and the boundary with information read from the datafile
         * Cuts the domain at x=0, creating a Left and a Right portion of the bulk
         */
        virtual void initRegions(getfem::mesh&) const override;

    };


    class GmshBuilder : public MeshBuilderStrategy {
        std::string M_meshFile; ///< The mesh filename

    public:
        
        GmshBuilder(const std::string& meshFile):
        M_meshFile(meshFile){}

        /**
         * @brief Builds the mesh reading information from the .msh file
         */
        void buildMesh(getfem::mesh&) const override;

        /**
         * @brief Builds the mesh regions with information provided in the .msh file
         */
        virtual void initRegions(getfem::mesh&) const override;
    };




}

#endif // _MESH_BUILDER_HPP_