#include "BCHandler.hpp"
#include "Utils.hpp"
#include <sstream>

namespace gf {

    BCHandler::BCHandler(const getfem::mesh& mesh, const Domain& domain)
    : M_mesh(mesh) {

        // Get the domain dimensions
        scalar_type Lx = domain.Lx;
        scalar_type Ly = domain.Ly; 
        scalar_type Lz = domain.Lz;

        // Construct the map regionID -> meshRegion
        getfem::outer_faces_of_mesh(mesh, M_border_faces);

        M_IdToRegion[1] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{-Lx/2.,0.,0.},
            base_node{0.,Ly,0.});
        M_IdToRegion[2] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{0.,Ly,0.},
            base_node{-Lx/2,Ly,Lz});
        M_IdToRegion[3] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{-Lx/2,Ly,Lz},
            base_node{0.,0.,Lz});
        M_IdToRegion[4] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{0.,0.,Lz},
            base_node{-Lx/2.,0.,0.});
        M_IdToRegion[5] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{0.,0.,0.},
            base_node{Lx/2.,Ly,0.});
        M_IdToRegion[6] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{Lx/2.,Ly,0.},
            base_node{0.,Ly,Lz});
        M_IdToRegion[7] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{0.,Ly,Lz},
            base_node{Lx/2.,0.,Lz});
        M_IdToRegion[8] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{Lx/2.,0.,Lz},
            base_node{0.,0.,0.});

        // alternatively:
        // M_IdToRegion[9] = getfem::select_faces_of_normal(mesh, border_faces,
        //    base_small_vector(-1.,0.,0.), 0.01);
        // M_IdToRegion[10] = getfem::select_faces_of_normal(mesh, border_faces,
        //    base_small_vector(1.,0.,0.), 0.01);
        M_IdToRegion[9] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{-Lx/2.,0.,0.},
            base_node{-Lx/2.,Ly,Lz});
        M_IdToRegion[10] = getfem::select_faces_in_box(
            mesh, M_border_faces,
            base_node{Lx/2.,0.,0},
            base_node{Lx/2.,Ly,Lz});


        // ALTERNATIVE:
        // for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
        //     assert(it.is_face());
        //     base_node x = /* get the centroid of the face*/

        //     base_small_vector un = mesh.normal_of_face_of_convex(it.cv(), it.f());
        //     un /= gmm::vect_norm2(un);

        //     if ((un[2] + 1.0) < 1.0E-7)
        //         M_IdToRegion[1].add(it.cv(), it.f());
        //     else if ((un[1] - 1.0 < 1.0E-7) && x[0] < 0) 
        //         M_IdToRegion[2].add(it.cv(),it.f());
        //     else if ((un[2] - 1.0 < 1.0E-7) && x[0] < 0)
        //         M_IdToRegion[3].add(it.cv(),it.f());
        //     else if ((un[1] + 1.0 < 1.0E-7) && x[0] < 0)
        //         M_IdToRegion[4].add(it.cv(),it.f());
        //     else if ((un[2] + 1.0) < 1.0E-7)
        //         M_IdToRegion[5].add(it.cv(), it.f());
        //     else if ((un[1] - 1.0 < 1.0E-7) && x[0] > 0) 
        //         M_IdToRegion[6].add(it.cv(),it.f());
        //     else if ((un[2] - 1.0 < 1.0E-7) && x[0] > 0)
        //         M_IdToRegion[7].add(it.cv(),it.f());
        //     else if ((un[1] + 1.0 < 1.0E-7) && x[0] > 0)
        //         M_IdToRegion[8].add(it.cv(),it.f());
        //     else if (un[0] + 1.0 < 1.0E-7)
        //         M_IdToRegion[9].add(it.cv(),it.f());
        //     else if (un[0] - 1.0 < 1.0E-7)
        //         M_IdToRegion[10].add(it.cv(),it.f());
        // }
    }

    void BCHandler::readBC(GetPot& datafile) {
        getfem::outer_faces_of_mesh(M_mesh, M_border_faces);
        read<BCType::Dirichlet>(datafile);
        read<BCType::Neumann>(datafile);
        read<BCType::Mixed>(datafile);
    }

    VectorFunctionType
    BCHandler::buildBCFunctionFromExpressions(const std::vector<std::string>& components)
    {
        // RMK: i assume that each call to buildBCFunctionFromExpressions resets the parser
        
        // Build a function for each component and return the combined VectorFunctionType
        VectorFunctionType parsedVectorFunction;
        std::vector<ScalarFunctionType> parsedFunctionsVec;

        for (size_type k {}; k < 3; ++k) {
            // Reset the expression of the parser for the new component
            M_parser.set_expression(components[k]);

            // Define a lambda function that binds to the parsed expression
            auto func = [this](base_node x, scalar_type t) -> scalar_type {
                // evaluate the expression for x[0], x[1], x[2], t and return the result as a base_small_vector
                std::array<double, 4> inputs = { x[0], x[1], x[2], t }; // assuming base_node has x[0], x[1], x[2]
                scalar_type result = M_parser(inputs);  // evaluating the function at the input values
                return result;  // return result as scalar_type
            };

            parsedFunctionsVec[k] = std::move(func);
        }

        // Combine and return the result
        return [parsedFunctionsVec](base_node node, scalar_type t) -> base_small_vector {
            base_small_vector result(3);
            for (size_type i {}; i < 3 ; ++i) {
                result[i] = parsedFunctionsVec[i](node, t);
            }
            return result;
        };

        return parsedVectorFunction;
    }


    template <BCType T>
    void BCHandler::read(GetPot& datafile){

        std::vector<size_type> regionsID;
        std::string regionsStr;
        if constexpr (T == BCType::Dirichlet) { // read the regionDisp list
            regionsStr = datafile("physics.regionDisp", "");
        }
        else if constexpr (T == BCType::Neumann) {
            regionsStr = datafile("physics.regionLoad", "");
        }
        else if constexpr (T == BCType::Mixed) {
            regionsStr = datafile("physics.regionMix", "");
        }

        std::istringstream istr(regionsStr);
        size_type val;
        while (istr >> val) regionsID.emplace_back(val);

        for (size_t i = 0; i < regionsID.size(); ++i) {
            std::ostringstream varname;

            if constexpr (T == BCType::Dirichlet) { // read the bdDisp list
                varname << "physics.bdDisp" << (i + 1);  // bdDisp1, bdDisp2, ...
                if (!datafile.search(varname.str().c_str()))
                    throw std::logic_error("Missing Dirichlet expression");
            }
            else if constexpr (T == BCType::Neumann) {
                varname << "physics.bdLoad" << (i + 1);  // bdLoad1, bdLoad2, ...
                if (!datafile.search(varname.str().c_str()))
                    throw std::logic_error("Missing Neumann expression");
            }
            else if constexpr (T == BCType::Mixed) {
                 /** !\todo */
            }

            std::string stringValue = datafile(varname.str().c_str(), "");

            std::vector<std::string> components = splitString(stringValue);

            // Parse using muParserInterface
            VectorFunctionType func = buildBCFunctionFromExpressions(components);

            if constexpr (T == BCType::Dirichlet) { // build BCDir and add to M_BCList
                // Build the BCDir object bc
                auto bc = std::make_unique<BCDir>(M_IdToRegion[regionsID[i]], regionsID[i], std::move(func), BCType::Dirichlet);
                // Add to map
                M_BCList[BCType::Dirichlet].emplace_back(std::move(bc));
            }
            else if constexpr (T == BCType::Neumann) {
                // Build the BCNeu object bc
                auto bc = std::make_unique<BCNeu>(M_IdToRegion[regionsID[i]], regionsID[i], std::move(func), BCType::Neumann);
                // Add to map
                M_BCList[BCType::Neumann].emplace_back(std::move(bc));
            }
            // else if constexpr (T == BCType::Mixed) { /** !\todo */
            //     // Build the BCMixed object bc
            //     auto bc = std::make_unique<BCMix>(M_IdToRegion[regionsID[i]], regionsID[i], std::move(func) /** !\todo */);
            //     // Add to map
            //     M_BCList[BCType::Mixed].emplace_back(std::move(bc));
            // }

        }

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



}