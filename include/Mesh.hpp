#ifndef _MESH_HPP_
#define _MESH_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "MeshBuilder.hpp"
#include "MeshRegion.hpp"
#include "GetPot"

namespace gf {

    class Mesh {
        const Params& M_params; ///< the Parameters
        std::unique_ptr<MeshBuilderStrategy> M_meshBuilder; ///< MeshBuilder
        getfem::mesh M_mesh; ///< The mesh
        RegionMapType M_regions; ///< The mesh regions (BulkLeft, BulkRight, Fault)
        BoundaryMapType M_bdRegions; /// Map (ID, BoundaryRegion)
        

    public:
        Mesh (const Params&);

        const getfem::mesh& get() const { return M_mesh; }

        const RegionMapType& getRegions() const { return M_regions; }

        const BoundaryMapType& getBdRegions() const { return M_bdRegions; }
        

    };

} // namespace gf

#endif // _MESH_HPP_