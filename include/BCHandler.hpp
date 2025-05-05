#ifndef _BC_HANDLER_HPP_
#define _BC_HANDLER_HPP_

#include "Core.hpp"
#include "GetPot"
#include "BC.hpp"

#include <unordered_map>

namespace gf {

    class BCHandler {

        using BCListType = std::unordered_map<BCType, std::vector<unique_ptr<BC>>>;

        const getfem::mesh& M_mesh; ///< The mesh object
        getfem::mesh_region M_border_faces; ///< Region containing all the border faces of the mesh
        std::unordered_map<size_type,getfem::mesh_region> M_IdToRegion;
        BCListType M_BCList;
        
        MuParserInterface::muParserXInterface<4> M_parser;

        VectorFunctionType buildBCFunctionFromExpressions(const std::vector<std::string>&);

        template <BCType T>
        void read(const GetPot&);

    public:
        BCHandler(const getfem::mesh& mesh, const Domain& domain);

        void readBC(const GetPot &);

    /*
        std::vector<base_node> getDirichletNodes() const;
        std::vector<base_node> getNeumannNodes() const;
        std::vector<base_node> getMixedNodes() const;
        std::vector<size_type> getDirichletRegions() const;
        std::vector<size_type> getNeumannRegions() const;
        std::vector<size_type> getMixedRegions() const;
    */
    }

} // namespace gf



#endif // _BC_HANDLER_HPP_