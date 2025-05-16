#ifndef _BC_HANDLER_HPP_
#define _BC_HANDLER_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "Mesh.hpp"
#include "GetPot"
#include "BC.hpp"
#include "muParserXInterface.hpp"

#include <unordered_map>

namespace gf {

    class BCHandler {

        using BCListType = std::unordered_map<BCType, std::vector<std::unique_ptr<BC>>>;
        using BCstringsType = std::unordered_map<BCType, std::vector<std::string>>;

        const getfem::mesh& M_mesh; ///< The mesh object
        BCListType M_BCList; ///< The BCobjects
        BCstringsType M_BCStrings; ///< The strings of the prescribed BC functions (just for output)
        
        muParserXInterface M_parser;

        // VectorFunctionType buildBCFunctionFromExpressions(const std::vector<std::string>&);

        template <BCType T>
        void read(const GetPot&, const BoundaryMapType&);

    public:
        BCHandler(const getfem::mesh&);

        void readBC(const GetPot &, const BoundaryMapType&);

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