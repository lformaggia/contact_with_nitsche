
#ifndef _MESH_BUILDER_HPP_
#define _MESH_BUILDER_HPP_

#include "Core.hpp"

namespace gf {

    class MeshRegion;

    /**
     * Builder class for building the mesh
     */
    class MeshBuilderStrategy{
    public:

        using RegionMapType = std::map<std::string, std::unique_ptr<MeshRegion>>;
        /**
         * @brief Null method to build the mesh, overridden in the hierarchy
         * They get a reference to the mesh object, modifying it
         */
        virtual void buildMesh(getfem::mesh&) const = 0;

        virtual void initRegions(const getfem::mesh&, RegionMapType&) const = 0;

        ~MeshBuilderStrategy() = default;
    };

    
    class GmshBuilder : public MeshBuilderStrategy {
        std::string M_meshFile;

    public:
        GmshBuilder(const std::string& meshFile):
        M_meshFile(meshFile){}

        void buildMesh(getfem::mesh&) const override;

        virtual void initRegions(const getfem::mesh&, RegionMapType&) const override;
    };

    class BuiltInBuilder : public MeshBuilderStrategy {
    private:
        GetPot M_datafile;

    public:
        BuiltInBuilder(const GetPot& datafile):
        M_datafile(datafile){}

        void buildMesh(getfem::mesh&) const override;

        virtual void initRegions(const getfem::mesh&, RegionMapType&) const override;

    };

}

#endif // _MESH_BUILDER_HPP_