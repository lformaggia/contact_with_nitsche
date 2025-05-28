#ifndef _MESH_HPP_
#define _MESH_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "MeshBuilder.hpp"
#include "GetPot"

namespace gf {

    class Mesh {

        const Params& M_params; ///< the Parameters
        std::unique_ptr<MeshBuilderStrategy> M_meshBuilder; ///< MeshBuilder
        getfem::mesh M_mesh; ///< The mesh

    public:
        /** Constructor*/
        Mesh (const Params&);

        /** @brief Return a const reference to the getfem::mesh object */
        const getfem::mesh& get() const { return M_mesh; }

        /** @brief Return the const region r */
        const getfem::mesh_region region(RegionType r) const { return M_mesh.region(r); }
        
        /** @brief Return the const boundary region i */
        const getfem::mesh_region bdRegion(size_type i) const { return M_mesh.region(i); }

        size_type dim() const { return M_mesh.dim(); }

    };

} // namespace gf

#endif // _MESH_HPP_