#include "BCHandler.hpp"
#include "Utils.hpp"
#include <sstream>
#include <algorithm>

bool DEBUGBC = false;

namespace gf {

    BCHandler::BCHandler(const getfem::mesh& m)
    : M_mesh(m) {
    }

    void BCHandler::readBC(const GetPot& gp) {
        read<BCType::Dirichlet>(gp);
        read<BCType::Neumann>(gp);
        read<BCType::Mixed>(gp);
    }


    template <BCType T>
    void BCHandler::read(const GetPot& datafile){

        std::string regionsStr;
        std::string normalDir;
        if constexpr (T == BCType::Dirichlet) // read the regionDisp list
            regionsStr = datafile("physics/regionDisp", "");
        else if constexpr (T == BCType::Neumann) 
            regionsStr = datafile("physics/regionLoad", "");
        else if constexpr (T == BCType::Mixed) 
            regionsStr = datafile("physics/regionDispNormal", "");

        std::vector<std::size_t> regionsID = gf::toVec(regionsStr);
        std::vector<std::size_t> regionsIDnormalDisp = gf::toVec(normalDir);

        M_BCStrings[T].reserve(regionsID.size()); // not really needed

        for (size_t i = 0; i < regionsID.size(); ++i) {
            std::ostringstream varname;

            if constexpr (T == BCType::Dirichlet) // read the bdDisp list
                varname << "physics/bdDisp" << (i + 1); // bdDisp1, bdDisp2, ...
            else if constexpr (T == BCType::Neumann)
                varname << "physics/bdLoad" << (i + 1); // bdLoad1, bdLoad2, ...
            else if constexpr (T == BCType::Mixed) {
                varname << "physics/bdDispN" << (i + 1);
            }

            std::string stringValue = datafile(varname.str().c_str(), "");

            if (stringValue.empty())
                throw std::runtime_error("String Values undetected!");
            
            /** @remark Alternatively to muParserX... 
             * std::vector<std::string> components = splitString(stringValue); 
             * VectorFunctionType func = buildBCFunctionFromExpressions(components);
            */

            // Use muparser (alternative to buildBCFunctionFromExpressions)
            M_parser.set_expression(stringValue);
            M_BCStrings[T].emplace_back(stringValue);

            /** \todo: check if n is needed (eventually change BC constructor) */
            getfem::mr_visitor it(M_mesh.region(regionsID[i]));
            base_small_vector n = M_mesh.mean_normal_of_face_of_convex(it.cv(), it.f());

            if constexpr (T == BCType::Dirichlet) { // build BCDir and add to M_BCList
                // Build the BCDir object bc
                auto bc = std::make_unique<BCDir>(M_mesh.region(regionsID[i]), regionsID[i], M_parser, T, n);
                // Add to map
                M_BCList[T].emplace_back(std::move(bc));
            }

            else if constexpr (T == BCType::Neumann) {
                // Build the BCNeu object bc
                auto bc = std::make_unique<BCNeu>(M_mesh.region(regionsID[i]), regionsID[i], M_parser, T, n);
                // Add to map
                M_BCList[T].emplace_back(std::move(bc));
            }

            else if constexpr (T == BCType::Mixed) { /** !\todo */
                // Build the BCMixed object bc
                auto bc = std::make_unique<BCMix>(M_mesh.region(regionsID[i]), regionsID[i], M_parser, T, n);
                // Add to map
                M_BCList[BCType::Mixed].emplace_back(std::move(bc));
            }

            if (DEBUGBC)
                if constexpr(T==BCType::Dirichlet)
                std::clog << "DEBUG: evaluating bdDisp" << i <<" at ((1,2,3),100): ("
                << M_BCList[BCType::Dirichlet][i]->eval({1,2,3},100)[0] << ", "
                << M_BCList[BCType::Dirichlet][i]->eval({1,2,3},100)[1] << ", "
                << M_BCList[BCType::Dirichlet][i]->eval({1,2,3},100)[2]
                << ")" << std::endl;

        }

    }

    const std::vector<std::unique_ptr<BC>> & 
    BCHandler::Neumann() const {
        static const std::vector<std::unique_ptr<BC>> empty;
        auto it = M_BCList.find(BCType::Neumann);
        return (it != M_BCList.end()) ? it->second : empty;
        // return M_BCList.at(BCType::Neumann);
    }

    const std::vector<std::unique_ptr<BC>> &
    BCHandler::Dirichlet() const {
        static const std::vector<std::unique_ptr<BC>> empty;
        auto it = M_BCList.find(BCType::Dirichlet);
        return (it != M_BCList.end()) ? it->second : empty;
        // return M_BCList.at(BCType::Dirichlet);
    }

    const std::vector<std::unique_ptr<BC>> & 
    BCHandler::Mixed() const {
        static const std::vector<std::unique_ptr<BC>> empty;
        auto it = M_BCList.find(BCType::Mixed);
        return (it != M_BCList.end()) ? it->second : empty;
    }

} // namespace gf