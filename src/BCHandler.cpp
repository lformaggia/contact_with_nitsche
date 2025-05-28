#include "BCHandler.hpp"
#include "Utils.hpp"
#include <sstream>
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
        if constexpr (T == BCType::Dirichlet) // read the regionDisp list
            regionsStr = datafile("physics/regionDisp", "");
        else if constexpr (T == BCType::Neumann) 
            regionsStr = datafile("physics/regionLoad", "");
        else if constexpr (T == BCType::Mixed)
            regionsStr = datafile("physics/regionMix", "");

        std::vector<std::size_t> regionsID = gf::toVec(regionsStr);
        std::cout << "Region string" << regionsStr << std::endl;
        std::cout << "ToVec result: ";
        for (auto el: regionsID) std::cout << el << " ";
        std::cout << std::endl;

        M_BCStrings[T].reserve(regionsID.size()); // not really needed

        for (size_t i = 0; i < regionsID.size(); ++i) {
            std::ostringstream varname;

            if constexpr (T == BCType::Dirichlet) // read the bdDisp list
                varname << "physics/bdDisp" << (i + 1); // bdDisp1, bdDisp2, ...
            else if constexpr (T == BCType::Neumann)
                varname << "physics/bdLoad" << (i + 1); // bdLoad1, bdLoad2, ...
            else if constexpr (T == BCType::Mixed) {
                 /** !\todo */
            }

            std::string stringValue = datafile(varname.str().c_str(), "");
            std::cout << varname.str() << std::endl;
            if (stringValue.empty())
                throw std::runtime_error("String Values undetected!");

            /** @remark Alternatively to muParserX... 
             * std::vector<std::string> components = splitString(stringValue); 
             * VectorFunctionType func = buildBCFunctionFromExpressions(components);
            */

            // Use muparser (alternative to buildBCFunctionFromExpressions)
            M_parser.set_expression(stringValue);
            M_BCStrings[T].emplace_back(stringValue);

            if constexpr (T == BCType::Dirichlet) { // build BCDir and add to M_BCList
                // Build the BCDir object bc
                auto bc = std::make_unique<BCDir>(M_mesh.region(regionsID[i]), regionsID[i], M_parser, T);
                // Add to map
                M_BCList[T].emplace_back(std::move(bc));
            }

            else if constexpr (T == BCType::Neumann) {
                // Build the BCNeu object bc
                auto bc = std::make_unique<BCNeu>(M_mesh.region(regionsID[i]), regionsID[i], M_parser, T);
                // Add to map
                M_BCList[T].emplace_back(std::move(bc));
            }

            // else if constexpr (T == BCType::Mixed) { /** !\todo */
            //     // Build the BCMixed object bc
            //     auto bc = std::make_unique<BCMix>(bds.at(regionsID[i]), regionsID[i], std::move(func) /** !\todo */);
            //     // Add to map
            //     M_BCList[BCType::Mixed].emplace_back(std::move(bc));
            // }

        }

        if (DEBUGBC)
            if constexpr(T==BCType::Dirichlet)
            std::clog << "DEBUG: evaluating bdDisp" <<" at ((1,2,3),100): ("
            << M_BCList[BCType::Dirichlet][1]->eval({1,2,3},100)[0] << ", "
            << M_BCList[BCType::Dirichlet][1]->eval({1,2,3},100)[1] << ", "
            << M_BCList[BCType::Dirichlet][1]->eval({1,2,3},100)[2]
            << ")" << std::endl;

        // THIS WILL GO INTO THE ASSEMBLY PHASE
        // gmm::resize(RM, mf.nb_dof(), mf.nb_dof());
        // gmm::resize(B, mf.nb_dof());
        // std::vector<double> F(mf.nb_basic_dof() * mf.get_qdim());

        // // Set Dirichlet values on boundary
        // for (auto it = mf.basic_dof_on_region(boundary).begin();
        //     it != mf.basic_dof_on_region(boundary).end(); ++it) {
        //     for (size_type q = 0; q < mf.get_qdim(); ++q)
        //         F[*it * mf.get_qdim() + q] = 0.0; // Zero Dirichlet BC
        // }

        // assembling_Dirichlet_condition(RM, B, mf, boundary, F);

    }

} // namespace gf