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
        /**
         * @brief Constructor that initializes the mesh.
         * @param p The Params object containing mesh parameters.
         * @note This constructor creates a MeshBuilderStrategy based on the mesh type specified in Params.
         *       It then builds the mesh and initializes regions.
         * @throws std::runtime_error if the mesh type is not recognized.
         * @throws std::runtime_error if the mesh cannot be built.
         */
        Mesh (const Params&);

        /** @brief Return a const reference to the getfem::mesh object */
        const getfem::mesh& get() const { return M_mesh; }

        /** @brief Return the const region r */
        const getfem::mesh_region region(RegionType r) const { return M_mesh.region(r); }
        
        /** @brief Return the mesh dimension (currently only 3D is supported) */
        size_type dim() const { return M_mesh.dim(); }

    };

} // namespace gf

#endif // _MESH_HPP_