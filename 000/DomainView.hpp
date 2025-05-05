#ifndef _DOMAIN_VIEW_HPP_
#define _DOMAIN_VIEW_HPP_

#include "Core.hpp"
#include <map>

namespace gf{

    class Domain;


    class DomainView {
    private:

    protected:
        getfem::mesh& M_mesh; ///< reference to the Domain object
        getfem::mesh_region M_region; ///< the mesh_region object from GetFEM

    public:

        /**
         * @brief Constructor for imported gmsh file
         * @param domain The Domain object to be referenced in the View
         */
        DomainView(getfem::mesh& mesh)
        : M_mesh(mesh) {}

        /**
         * @brief Constructor for imported gmsh file
         * @param domain The Domain object to be referenced in the View
         * @param region The physical region
         */
        DomainView(getfem::mesh& mesh, const getfem::mesh_region& region)
        : M_mesh(mesh), M_region(region) {}

        /**
         * @brief Constructor for built-in mesh
         * @param domain The Domain object to be referenced in the View
         * @param convexesList TList of elements in the region, with different behviours for the derived
         */
        DomainView(getfem::mesh& mesh, const dal::bit_vector& convexesList)
        : M_mesh(mesh), M_region(convexesList) {}
        
        virtual ~DomainView() = default;

        virtual std::string name() const = 0;

        const getfem::mesh_region& region() const { return M_region; }
        // const getfem::mesh& mesh() const; // via M_domain
    };


    class BulkView : public DomainView {
        SideType M_side;
    public:
        /**
         * @brief Constructs the BulkView passing the list of convexes belonging to the BulkRegion
         * This overload is used for a built-in mesh
         * It delegates to the base constructor
         */
        BulkView(getfem::mesh& mesh, const dal::bit_vector& convexesList, SideType s)
        : DomainView(mesh, convexesList), M_side(s){}

        /**
         * @brief Construct the BulkView importing from gmsh file
         * !\todo
         */
        BulkView(getfem::mesh& mesh, const getfem::mesh_region& region, SideType s)
        : DomainView(mesh,region), M_side(s){}

        /**
         * @brief Return the name of the region
         */
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


    class FaultView : public DomainView {
    public:
        using faceToConvexesMap = std::map<size_type, std::pair<size_type,size_type>>;
        /**
         * @brief Constructs the FaultView passing the list of convexes having points at x = 0 (interface)
         * This overload is used for a built-in mesh
         */
        FaultView(getfem::mesh& mesh, const dal::bit_vector&);

        /**
         * @brief Constructs the FaultView importing from gmsh file
         * !\todo: build the faceToConvexes map
         */
        FaultView(getfem::mesh& mesh, const getfem::mesh_region& region)
        : DomainView(mesh,region) {}

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

}

#endif // _DOMAIN_VIEW_HPP_
