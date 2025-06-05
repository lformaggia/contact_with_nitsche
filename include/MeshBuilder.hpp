
#ifndef _MESH_BUILDER_HPP_
#define _MESH_BUILDER_HPP_

#include "Core.hpp"
#include "Params.hpp"

namespace gf {

    class MeshRegion;

    /**
     * @brief Builder class for building the mesh
     * In the Mesh class constructor, the MeshBuilderStrategy::construct() method
     * is called. Depending on the command line argument passed, it either builds
     * the mesh internally or imports it from gmsh
     */
    class MeshBuilderStrategy{
    public:

        /**
         * @brief Constructor taking the domain parameters
         */
        MeshBuilderStrategy(const Domain& d)
        : M_domain(d) {}

        /**
         * @brief Construct the mesh by building it and initializing regions.
         * This is the entry point for all mesh creation logic.
         */
        void construct(getfem::mesh&) const;

        virtual ~MeshBuilderStrategy() = default;

    protected:

        Domain M_domain;

        /**
         * @brief Null method to build the mesh, overridden in the hierarchy
         * They get a reference to the mesh object, modifying it
         */
        virtual void buildMesh(getfem::mesh&) const = 0;

        /**
         * @brief Initializes the region
         * The method does nothing by default, and it's overridden by the
         * built-in builder only
         */
        virtual void initRegions(getfem::mesh&) const = 0;

    };


    class BuiltInBuilder : public MeshBuilderStrategy {

    public:
        BuiltInBuilder(const Domain& d):
        MeshBuilderStrategy(d) {}

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

    public:
        
        GmshBuilder(const Domain& d):
        MeshBuilderStrategy(d) {}

        /**
         * @brief Builds the mesh reading information from the .msh file
         */
        void buildMesh(getfem::mesh&) const override;

        /**
         * @brief Post-process the fault regions automatically imported
         */
        virtual void initRegions(getfem::mesh&) const override;
    
    private:
        
        void generateMeshFile() const;

    };




}

#endif // _MESH_BUILDER_HPP_