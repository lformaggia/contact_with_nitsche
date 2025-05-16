#ifndef _MESH_REGION_HPP_
#define _MESH_REGION_HPP_

#include "Core.hpp"
#include <map>

namespace gf{

    class MeshRegion {
    private:

    protected:
        const getfem::mesh& M_mesh; ///< reference to the Domain object
        getfem::mesh_region M_region; ///< the mesh_region object from GetFEM

    public:

        /**
         * @brief Constructor for imported gmsh file
         * @param mesh The mesh object
         */
        MeshRegion(const getfem::mesh& mesh)
        : M_mesh(mesh) {}

        /**
         * @brief Constructor for imported gmsh file
         * @param mesh The mesh object
         * @param region The physical region
         */
        MeshRegion(const getfem::mesh& mesh, const getfem::mesh_region& region)
        : M_mesh(mesh), M_region(region) {}

        /**
         * @brief Constructor for built-in mesh
         * @param domain The Domain object to be referenced in the View
         * @param convexesList TList of elements in the region, with different behviours for the derived
         */
        MeshRegion(const getfem::mesh& mesh, const dal::bit_vector& convexesList)
        : M_mesh(mesh), M_region(convexesList) {}
        
        virtual ~MeshRegion() = default;

        /**
         * @brief Return the name of the region
         */
        virtual std::string name() const = 0;

        /**
         * @brief Return a const& to the region
         */
        const getfem::mesh_region& region() const { return M_region; }

        /**
         * @brief Return the convex indeices of the region
         */
        const dal::bit_vector& index() const { return M_region.index(); }
    };


    class Bulk : public MeshRegion {
        SideType M_side;
    public:
        /**
         * @brief Constructs the Bulk passing the list of convexes belonging to the BulkRegion
         * This overload is used for a built-in mesh
         * It delegates to the base constructor
         */
        Bulk(const getfem::mesh& mesh, const dal::bit_vector& convexesList, SideType s)
        : MeshRegion(mesh, convexesList), M_side(s){}

        /**
         * @brief Construct the Bulk importing from gmsh file
         * !\todo
         */
        Bulk(const getfem::mesh& mesh, const getfem::mesh_region& region, SideType s)
        : MeshRegion(mesh,region), M_side(s){}

        std::string name() const override {
            if (M_side == SideType::LEFT) return "BulkLeft";
            return "BulkRight";
        }

        /**
         * @brief Return the list of convexes with a face lying on the interface
         */
        const dal::bit_vector interfConvex() const;

        /*...*/


    };


    class Fault : public MeshRegion {
    public:
        using faceToConvexesMap = std::map<size_type, std::pair<size_type,size_type>>;
        /**
         * @brief Constructs the Fault passing the list of convexes having points at x = 0 (interface)
         * This overload is used for a built-in mesh
         */
        Fault(const getfem::mesh& mesh, const dal::bit_vector&);

        /**
         * @brief Constructs the Fault importing from gmsh file
         * !\todo: build the faceToConvexes map
         */
        Fault(const getfem::mesh& mesh, const getfem::mesh_region& region)
        : MeshRegion(mesh,region) {}

        std::string name() const override { return "Fault"; }

        const std::vector<base_node> getInterfaceNodes() const;

        faceToConvexesMap faceToConvexes() const { return convexesSharingFace; }

    private:
        base_small_vector normal;
        base_small_vector tangent;
        faceToConvexesMap convexesSharingFace;

        /**
         * !\todo: add the following methods/members:
         * 
         */

    };

    class Boundary : public MeshRegion {
        size_type M_ID;
    public:
        /**
         * @brief Build the map (ID, BoundaryRegion) needed to assembly the boundary conditions
         */
        Boundary(const getfem::mesh& mesh, const getfem::mesh_region& region, size_type id)
        : MeshRegion(mesh,region), M_ID(id){}

        std::string name() const override {
            std::stringstream ss;
            ss << "Boundary" << M_ID;
            return ss.str();
        }

    };

}

#endif // _MESH_REGION_HPP_
