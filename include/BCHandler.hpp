#ifndef _BC_HANDLER_HPP_
#define _BC_HANDLER_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "GetPot"
#include "BC.hpp"
#include "muParserXInterface.hpp"

#include <unordered_map>

namespace gf {

    class BCHandler {

        using BCListType = std::unordered_map<BCType, std::vector<std::unique_ptr<BC>>>;

        const getfem::mesh& M_mesh; ///< The mesh object
        getfem::mesh_region M_border_faces; ///< Region containing all the border faces of the mesh
        std::unordered_map<size_type,getfem::mesh_region> M_IdToRegion;
        BCListType M_BCList;
        
        MuParserInterface2::muParserXInterface M_parser;

        VectorFunctionType buildBCFunctionFromExpressions(const std::vector<std::string>&);

        template <BCType T>
        void read(GetPot&);

    public:
        BCHandler() = default;
        BCHandler(const getfem::mesh& mesh, const Domain& domain);

        void readBC(GetPot &);

    /*
        std::vector<base_node> getDirichletNodes() const;
        std::vector<base_node> getNeumannNodes() const;
        std::vector<base_node> getMixedNodes() const;
        std::vector<size_type> getDirichletRegions() const;
        std::vector<size_type> getNeumannRegions() const;
        std::vector<size_type> getMixedRegions() const;
    */
    };

} // namespace gf



#endif // _BC_HANDLER_HPP_