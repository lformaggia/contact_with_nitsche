#include "BCHandler.hpp"
#include "Utils.hpp"
#include <sstream>
bool DEBUGBC = false;

namespace gf {

    BCHandler::BCHandler(const getfem::mesh& m)
    : M_mesh(m) {
    }

    void BCHandler::readBC(const GetPot& gp, const BoundaryMapType& bds) {
        read<BCType::Dirichlet>(gp, bds);
        read<BCType::Neumann>(gp,bds);
        read<BCType::Mixed>(gp,bds);
    }

    // VectorFunctionType
    // BCHandler::buildBCFunctionFromExpressions(const std::vector<std::string>& components)
    // {
    //     // RMK: i assume that each call to buildBCFunctionFromExpressions resets the parser
        
    //     // Build a function for each component and return the combined VectorFunctionType
    //     VectorFunctionType parsedVectorFunction;
    //     std::vector<ScalarFunctionType> parsedFunctionsVec;

    //     for (size_type k {}; k < 3; ++k) {
    //         // Reset the expression of the parser for the new component
    //         M_parser.set_expression(components[k]);

    //         // Define a lambda function that binds to the parsed expression
    //         auto func = [this](base_node x, scalar_type t) -> scalar_type {
    //             // evaluate the expression for x[0], x[1], x[2], t and return the result as a base_small_vector
    //             std::array<double, 4> inputs = { x[0], x[1], x[2], t }; // assuming base_node has x[0], x[1], x[2]
    //             scalar_type result = M_parser(inputs);  // evaluating the function at the input values
    //             return result;  // return result as scalar_type
    //         };

    //         parsedFunctionsVec[k] = std::move(func);

    //     }

    //     // Combine and return the result
    //     return [parsedFunctionsVec](base_node node, scalar_type t) -> base_small_vector {
    //         base_small_vector result(3);
    //         for (size_type i {}; i < 3 ; ++i) {
    //             result[i] = parsedFunctionsVec[i](node, t);
    //         }
    //         return result;
    //     };

    //     return parsedVectorFunction;
    // }


    template <BCType T>
    void BCHandler::read(const GetPot& datafile, const BoundaryMapType& bds){

        std::string regionsStr;
        if constexpr (T == BCType::Dirichlet) // read the regionDisp list
            regionsStr = datafile("physics/regionDisp", "");
        else if constexpr (T == BCType::Neumann) 
            regionsStr = datafile("physics/regionLoad", "");
        else if constexpr (T == BCType::Mixed)
            regionsStr = datafile("physics/regionMix", "");

        std::vector<std::size_t> regionsID = gf::toVec(regionsStr);

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
                auto bc = std::make_unique<BCDir>(*(bds.at(regionsID[i])), regionsID[i], M_parser, T);
                // Add to map
                M_BCList[T].emplace_back(std::move(bc));
            }

            else if constexpr (T == BCType::Neumann) {
                // Build the BCNeu object bc
                auto bc = std::make_unique<BCNeu>(*(bds.at(regionsID[i])), regionsID[i], M_parser, T);
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