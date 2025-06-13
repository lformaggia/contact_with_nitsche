#include "ContactProblem.hpp"

namespace gf{

    bool contact_flag = false;

    ContactProblem::ContactProblem(const Mesh& m, const Params& p)
    : M_mesh(m),
    M_params(p),
    M_BC(m.get()),
    M_FEM(m.get()),
    M_integrationMethod(m.get()),
    M_imData(M_integrationMethod)
    {
        if (M_params.contact.method == "nitsche")
            M_contactEnforcement = std::make_unique<NitscheContactEnforcement>(
                M_params.contact.theta, M_params.contact.gamma0);
        else if (M_params.contact.method == "penalty")
            M_contactEnforcement = std::make_unique<PenaltyContactEnforcement>(
                M_params.contact.epsilon);
    }

    void
    ContactProblem::init()
    {

        M_BC.readBC(M_params.datafile);

        M_FEM.setMeshFem(M_params.numerics, M_mesh.get());

        getfem::pintegration_method ppi = getfem::int_method_descriptor(M_params.numerics.integration);
        size_type N = M_params.domain.dim;
        M_integrationMethod.set_integration_method(M_mesh.get().convex_index(), ppi);
    
    }


    void
    ContactProblem::assemble()
    {
        std::cout << "Preparing the assembly phase:\n";
        
        gmm::set_traces_level(1);
        
        dim_type dim = M_mesh.get().dim();
        size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();

        // Main unknown of the problem (displacement)
        std::cout << "  Defining variables...";
        
        M_model.add_fem_variable("uL", M_FEM.mf_u1());
        M_model.add_fem_variable("uR", M_FEM.mf_u2());
        
        std::cout << "done.\n";


        // Add scalar data to the model
        std::cout << "  Initializing scalar data...";
        
        M_model.add_initialized_scalar_data("lambda", M_params.physics.M_lambda);
        M_model.add_initialized_scalar_data("mu", M_params.physics.M_mu);
        M_model.add_initialized_scalar_data("mu_fric", M_params.physics.M_mu_friction);
        
        std::cout << "done.\n";


        // Define some useful macro for Nitsche contact integrals
        std::cout << "  Defining macros...";

        M_model.add_initialized_scalar_data("eps", 1.e-20);
        M_model.add_macro("n", "Normal"); // use normal on contact face
        M_model.add_macro("eps_dir", "1e-4"); // to avoid numerical degeneracies
        M_model.add_macro("a", "[0, 1, 0]"); /** \todo Modify a to be more robust */

        // Project a onto plane orthogonal to n
        M_model.add_macro("t1_raw", "a - (a . n) * n");
        M_model.add_macro("t1", "Normalized(t1_raw)");

        // Cross to get the second tangent
        M_model.add_macro("t2", "Cross_product(n, t1)");

        // Define macros for jump
        M_model.add_macro("u_jump", "(uL - Interpolate(uR,neighbor_element))");
        M_model.add_macro("v_jump", "(Test_uL - Interpolate(Test_uR,neighbor_element))");
        M_model.add_macro("un_jump", "u_jump . n"); // normal trial jump
        M_model.add_macro("ut1_jump", "u_jump . t1");
        M_model.add_macro("ut2_jump", "u_jump . t2");
        M_model.add_macro("vn_jump", "v_jump . n"); // normal test jump
        M_model.add_macro("vt1_jump", "v_jump . t1");
        M_model.add_macro("vt2_jump", "v_jump . t2");

        
        M_model.add_macro("stressL","(lambda*Trace(Grad_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_uL))");
        M_model.add_macro("stressR","(lambda*Trace(Grad_uR)*Id(qdim(uR)) + 2*mu*Sym(Grad_uR))");
        M_model.add_macro("stressL_voigt", "[stressL(1,1), stressL(2,2), stressL(3,3), 2*stressL(2,3), 2*stressL(1,3), 2*stressL(1,2)]");
        M_model.add_macro("stressR_voigt", "[stressR(1,1), stressR(2,2), stressR(3,3), 2*stressR(2,3), 2*stressR(1,3), 2*stressR(1,2)]");

        
        M_model.add_macro("traction_u", "((lambda*Trace(Grad_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_uL)) * n)");
        M_model.add_macro("traction_v", "((lambda*Trace(Grad_Test_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_Test_uL)) * n)");
        M_model.add_macro("sig_u_nL", "traction_u . n");
        M_model.add_macro("sig_u_t1", "traction_u . t1");
        M_model.add_macro("sig_u_t2", "traction_u . t2");
        M_model.add_macro("sig_v_nL", "traction_v . n");
        M_model.add_macro("sig_v_t1", "traction_v . t1");
        M_model.add_macro("sig_v_t2", "traction_v . t2");

        // M_model.add_initialized_scalar_data("gammaN", M_params.contact.gamma0);
        // M_model.add_initialized_scalar_data("theta", M_params.contact.theta);  // symmetric variant
        
        // // Normal gap and stress
        // M_model.add_macro("Pn_u", "(gammaN * un_jump - sig_u_nL)");
        // M_model.add_macro("Pt1_u", "(gammaN * ut1_jump - sig_u_t1)");
        // M_model.add_macro("Pt2_u", "(gammaN * ut2_jump - sig_u_t2)");
        // M_model.add_macro("Pn_v_theta", "(gammaN * vn_jump - theta*sig_v_nL)");
        // M_model.add_macro("Pt1_v_theta", "(gammaN * vt1_jump - theta*sig_v_t1)");
        // M_model.add_macro("Pt2_v_theta", "(gammaN * vt2_jump - theta*sig_v_t2)");
        
        // // Friction threshold
        // M_model.add_macro("Sh", "mu_fric * pos_part(Pn_u)");
        // M_model.add_macro("norm_Pt", "sqrt(Pt1_u*Pt1_u + Pt2_u*Pt2_u)");
        // M_model.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
        // M_model.add_macro("proj_Pt2_u", "Pt2_u * min(1, Sh / (norm_Pt + eps))");

        // std::cout << "done.\n";

        // Enforcement of contact conditions
        M_contactEnforcement->enforce(M_model, M_integrationMethod);

        // Add isotropic elasticity bricks
        std::cout << "  Adding elasticity bricks...";
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uL", "lambda", "mu", RegionType::BulkLeft);
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uR", "lambda", "mu", RegionType::BulkRight);
        std::cout << "done.\n";

        // // Add linear stress brick
        // std::cout << "  Adding linear stress brick...";

        // getfem::add_linear_term(
        //     M_model,
        //     M_integrationMethod,
        //     "- theta/gammaN * sig_u_nL * sig_v_nL", /** expression */
        //     Fault, /** region */
        //     false, /** symmetric */
        //     false, /** coercive */
        //     "linear_stress",
        //     false /** check */
        // );

        // std::cout << "done.\n";


        // // Add KKT condition brick
        // std::cout << "  Adding KKT condition brick...";

        // getfem::add_nonlinear_term(
        //     M_model,
        //     M_integrationMethod,
        //     "1/gammaN * pos_part(Pn_u) * Pn_v_theta",
        //     Fault,
        //     false,
        //     false,
        //     "KKTbrick"
        // );

        // std::cout << "done.\n";


        // // Add Coulomb condition brick
        // std::cout << "  Adding Coulomb friction brick...";
        
        // getfem::add_nonlinear_term(
        //     M_model,
        //     M_integrationMethod,
        //     "(1/gammaN) * (proj_Pt1_u * Pt1_v_theta + proj_Pt2_u * Pt2_v_theta)",
        //     Fault,
        //     false,
        //     false,
        //     "CoulombBrick"
        // );
        
        // std::cout << "done.\n";


        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "eps * uL . Test_uL + eps * uR . Test_uR",
            BulkLeft,
            false,
            false,
            "stabilizer1"
        );
        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "eps * uL . Test_uL + eps * uR . Test_uR",
            BulkRight,
            false,
            false,
            "stabilizer2"
        );

        
        // Volumic source term (gravity)
        std::cout << "  Adding volumic source term brick...";

        plain_vector G(M_FEM.mf_rhs().nb_dof()*dim);

        for (size_type i = 0; i < nb_dof_rhs; ++i)
            gmm::copy(M_params.physics.M_gravity, gmm::sub_vector(G, gmm::sub_interval(i*dim, dim)));

        M_model.add_initialized_fem_data("VolumicData", M_FEM.mf_rhs(), G);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uL", "VolumicData", BulkLeft);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uR", "VolumicData", BulkRight);
        
        std::cout << "done.\n";


        // Neumann conditions
        std::cout << "  Adding Neumann condition bricks...";
        
        const auto & NeumannBCs = M_BC.Neumann();
        plain_vector F(nb_dof_rhs*dim);
        for (const auto& bc: NeumannBCs){
            const auto& rg = bc->getRegion();
            auto& f = bc->f();
            // just for debug
            auto f0 = [&f](const base_node&x){ return f(x, 0); };
            getfem::interpolation_function(M_FEM.mf_rhs(),F,f0, rg);
            M_model.add_initialized_fem_data(bc->name(), M_FEM.mf_rhs(),F);

            if (bc->isLeftBoundary())
                getfem::add_source_term_brick(M_model, M_integrationMethod,"uL",bc->name(),bc->ID());
            else
                getfem::add_source_term_brick(M_model, M_integrationMethod,"uR",bc->name(),bc->ID());
        }
        std::cout << "done.\n";

        // Dirichlet conditions
        std::cout << "  Adding Dirichlet condition bricks...";
        const auto & DirichletBCs = M_BC.Dirichlet();
        plain_vector D(nb_dof_rhs*dim);

        for (const auto& bc: DirichletBCs){
            const auto& rg = bc->getRegion();
            auto& f = bc->f();
            // take time 0
            auto f0 = [&f](const base_node&x){ return f(x, 0); };

            /// @DEBUG:            
            getfem::interpolation_function(M_FEM.mf_rhs(),D,f0, rg);
            M_model.add_initialized_fem_data(bc->name(), M_FEM.mf_rhs(),D);

            if (bc->isLeftBoundary())
                getfem::add_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uL", M_FEM.mf_u1(), bc->ID(), 
                     bc->name());
            else
                getfem::add_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uR", M_FEM.mf_u2(), bc->ID(), 
                     bc->name());
        }

        std::cout << "done.\n";


        // Mixed conditions (normal Dirichlet)
        std::cout << "  Adding normal Dirichlet condition bricks...";

        const auto & MixedBCs = M_BC.Mixed();
        plain_vector M(nb_dof_rhs);

        for (const auto& bc: MixedBCs){
            const auto& rg = bc->getRegion();
            auto& f = bc->f();

            auto f0 = [&f](const base_node&x){ return f(x, 0); };
            getfem::interpolation_function(M_FEM.mf_rhs(), M, f0, rg);

            M_model.add_initialized_fem_data(bc->name(), M_FEM.mf_rhs(),M);

            if (bc->isLeftBoundary())
                // getfem::add_normal_Dirichlet_condition_with_penalization
                //     (M_model, M_integrationMethod, "uL", 1e-5, bc->ID(), bc->name());
                getfem::add_normal_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uL", M_FEM.mf_rhs(), bc->ID(), bc->name());
            else
                // getfem::add_normal_Dirichlet_condition_with_penalization
                //     (M_model, M_integrationMethod, "uR", 1e-5, bc->ID(), bc->name());
                getfem::add_normal_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uR", M_FEM.mf_rhs(), bc->ID(), bc->name());  
        }

        std::cout << "done.\n";

    }

    void
    ContactProblem::solve() {

        std::cout << "Solving the problem..." << std::endl;

        const auto & NeumannBCs = M_BC.Neumann();
        gmm::set_traces_level(1);
        
        dim_type dim = M_mesh.get().dim();
        size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();
        plain_vector F(nb_dof_rhs*dim);
        
        plain_vector VM(M_FEM.mf_stress().nb_dof());

        // Time loop
        size_type n_timesteps = static_cast<size_type>(M_params.time.tend - M_params.time.t0)/M_params.time.dt;
        
        scalar_type t {};

        for (size_type i{}; i < n_timesteps; ++i)
        {
            std::cout << "\nt = " << t << std::endl;
            // Neumann conditions
            for (const auto& bc: NeumannBCs){
                const auto& rg = bc->getRegion();
                auto& f = bc->f();
                // just for debug
                auto ft = [&f,t](const base_node&x){ return f(x, t); };
                getfem::interpolation_function(M_FEM.mf_rhs(),F,ft, rg);
                gmm::copy(F, M_model.set_real_variable(bc->name()));
            }

            // Solve the problem
            gmm::iteration iter(M_params.it.atol, 1, M_params.it.maxIt);
            getfem::standard_solve(M_model,iter);
            if (iter.converged()) {
                std::cout << "  Iteration converged after " << iter.get_iteration() << " iterations." << std::endl;
            } else {
                std::cerr << " Warning: Iteration did not converge after " << iter.get_iteration() << " iterations." << std::endl;
            }

            // Export results
            exportVtk(i);

            // Advance in time
            t += M_params.time.dt;
        }

    }

    
    void ContactProblem::exportVtk(size_type i) {

        dim_type dim = M_mesh.get().dim();

        // Zero out the noised components
        plain_vector uL_backup = M_model.real_variable("uL");
        plain_vector uR_backup = M_model.real_variable("uR");

        plain_vector &uL = M_model.set_real_variable("uL");
        plain_vector &uR = M_model.set_real_variable("uR");

        auto dofs_uL = M_FEM.mf_u1().dof_on_region(M_mesh.region(BulkLeft));
        auto dofs_uR = M_FEM.mf_u2().dof_on_region(M_mesh.region(BulkRight));
        auto dofs_uF = dofs_uL & dofs_uR; // shared dofs on the fault
        auto dofs_uL_filtered = dofs_uL;
        auto dofs_uR_filtered = dofs_uR;
        dofs_uL_filtered.setminus(dofs_uF);
        dofs_uR_filtered.setminus(dofs_uF);

        for (dal::bv_visitor i(dofs_uL_filtered); !i.finished(); ++i)
            uR[i] = 0.0;
        for (dal::bv_visitor j(dofs_uR_filtered); !j.finished(); ++j)
            uL[j] = 0.0;

        M_FEM.mf_stress().set_qdim(3,3); // Voigt: [Sxx, Syy, Szz, Syz, Sxz, Sxy]
        
        // Post-process stress
        std::cout << "Computing stress...";

        plain_vector VL(M_FEM.mf_stress().nb_dof());
        plain_vector VR(M_FEM.mf_stress().nb_dof());
        
        getfem::ga_interpolation_Lagrange_fem(M_model, "stressL", M_FEM.mf_stress(), VL);
        getfem::ga_interpolation_Lagrange_fem(M_model, "stressR", M_FEM.mf_stress(), VR);

        auto dofs_sL = M_FEM.mf_stress().dof_on_region(M_mesh.region(BulkLeft));
        auto dofs_sR = M_FEM.mf_stress().dof_on_region(M_mesh.region(BulkRight));
        auto dofs_sF = dofs_sL & dofs_sR; // shared dofs on the fault
        auto dofs_sL_filtered = dofs_sL;
        auto dofs_sR_filtered = dofs_sR;
        dofs_sL_filtered.setminus(dofs_sF);
        dofs_sR_filtered.setminus(dofs_sF);

        for (dal::bv_visitor i(dofs_sL_filtered); !i.finished(); ++i)
            VR[i] = 0.0;
        for (dal::bv_visitor j(dofs_sR_filtered); !j.finished(); ++j)
            VL[j] = 0.0;

        std::cout << "done." << std::endl;

        // Custom export to VTK
        std::cout << "Exporting results to result_" << i << ".vtk...";

        scalar_type offset = 1e-5; // offset for visualization

        // const plain_vector& uL = M_model.real_variable("uL");
        // const plain_vector& uR = M_model.real_variable("uR");

        size_type nb_dof_L = dofs_uL.card();
        size_type nb_dof_R = dofs_uR.card();
        size_type nb_dof_F = dofs_uF.card();
        size_type nb_dof = M_FEM.mf_u1().nb_dof();
        size_type nb_dof_tot = nb_dof_L + nb_dof_R; // == nb_dof + nb_dof_F

        plain_vector U;
        U.reserve(nb_dof_tot);
        std::vector<base_small_vector> U_vector;
        U_vector.reserve(U.size() / dim);

        std::vector<base_node> dof_coords;
        dof_coords.reserve(nb_dof_tot / dim);  // each node has 3 dofs

        // Map original dof index to new index in U
        std::map<size_type, size_type> dof_map;
        std::map<size_type, size_type> dof_to_point;

        // Insert uL values
        for (dal::bv_visitor ii(dofs_uL); !ii.finished(); ++ii) {
            size_type dof = ii;
            if (dof % dim == 0) {  // only once per node
                base_node pt = M_FEM.mf_u1().point_of_basic_dof(dof);
                pt[0] -= offset;
                dof_coords.push_back(pt);
                size_type pt_idx = dof_coords.size() - 1;
                dof_to_point[dof] = pt_idx;
            }
            U.push_back(uL[dof]);
        }

        // Insert uR values
        for (dal::bv_visitor jj(dofs_uR); !jj.finished(); ++jj) {
            size_type dof = jj;
            if (dof % dim == 0) {  // only once per node
                base_node pt = M_FEM.mf_u2().point_of_basic_dof(dof);
                pt[0] += offset;
                dof_coords.push_back(pt);
                size_type pt_idx = dof_coords.size() - 1;
                dof_to_point[dof] = pt_idx;
            }
            U.push_back(uR[dof]);
        }

        for (size_type i = 0; i < U.size(); i += dim) {
            base_small_vector vec(dim);
            for (size_type d = 0; d < dim; ++d)
                vec[d] = U[i + d];
            U_vector.emplace_back(vec);
        }

        assert(U_vector.size() == dof_coords.size());

        // Write stress in Voigt notation
        std::vector<std::array<double, 6>> stress_voigt(dof_coords.size());

        auto insert_stress_voigt = [&](const plain_vector &V, const dal::bit_vector &dofs, double offset_sign) {
            for (dal::bv_visitor ii(dofs); !ii.finished(); ++ii) {
                size_type dof = ii;
                if (dof % 9 == 0) {
                    base_node pt = M_FEM.mf_stress().point_of_basic_dof(dof);
                    pt[0] += offset_sign * offset;

                    for (size_type i = 0; i < dof_coords.size(); ++i) {
                        if (gmm::vect_dist2(pt, dof_coords[i]) < 1e-10) {
                            stress_voigt[i] = {
                                V[dof + 0],                          // Sxx
                                V[dof + 4],                          // Syy
                                V[dof + 8],                          // Szz
                                0.5 * (V[dof + 5] + V[dof + 7]),     // Syz
                                0.5 * (V[dof + 2] + V[dof + 6]),     // Sxz
                                0.5 * (V[dof + 1] + V[dof + 3])      // Sxy
                            };
                            break;
                        }
                    }
                }
            }
        };

        insert_stress_voigt(VL, dofs_sL, -1.0);
        insert_stress_voigt(VR, dofs_sR, +1.0);


        // Build connectivity
        std::vector<std::vector<size_type>> cells;
        std::vector<int> cell_types;

        auto reorderNodesForVTK = [](std::vector<size_type>& nodes, bool isHex) {
            if (isHex){
                std::swap(nodes[2], nodes[3]);
                std::swap(nodes[6], nodes[7]);
                return nodes;
            }
            // For tetrahedra, no reordering needed
            return nodes;
        };


        for (getfem::mr_visitor i(M_mesh.region(BulkLeft)); !i.finished(); ++i) {
            const bgeot::pconvex_structure cvs = M_mesh.get().structure_of_convex(i.cv());
            size_type nb_pts = cvs->nb_points();
            std::vector<size_type> pt_indices;

            for (size_type ipt = 0; ipt < nb_pts; ++ipt) {
                base_node pt = M_mesh.get().points_of_convex(i.cv())[ipt];
                // std::cout << "point " << ipt << " in cv " << i.cv() << ": "
                //     << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
                pt[0] -= offset;

                // Find corresponding index in dof_coords
                bool found = false;
                for (size_type i = 0; i < dof_coords.size(); ++i) {
                    if (gmm::vect_dist2(pt, dof_coords[i]) < 1e-10) {
                        pt_indices.push_back(i);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cerr << "Warning: point not found in coordinate list.\n";
                }
            }

            if (pt_indices.size() == nb_pts) {
                auto cvx_type = M_params.domain.meshType == "GT_PK(3,1)"? 10 : 12;
                cells.push_back(reorderNodesForVTK(pt_indices, cvx_type == 12));
                cell_types.push_back(cvx_type);
            }
        }

        for (getfem::mr_visitor i(M_mesh.region(BulkRight)); !i.finished(); ++i) {
            const bgeot::pconvex_structure cvs = M_mesh.get().structure_of_convex(i.cv());
            size_type nb_pts = cvs->nb_points();
            std::vector<size_type> pt_indices;

            for (size_type ipt = 0; ipt < nb_pts; ++ipt) {
                base_node pt = M_mesh.get().points_of_convex(i.cv())[ipt];
                // std::cout << "point " << ipt << " in cv " << i.cv() << ": "
                //     << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
                pt[0] += offset;

                // Find corresponding index in dof_coords
                bool found = false;
                for (size_type i = 0; i < dof_coords.size(); ++i) {
                    if (gmm::vect_dist2(pt, dof_coords[i]) < 1e-10) {
                        pt_indices.push_back(i);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cerr << "Warning: point not found in coordinate list.\n";
                }
            }

            if (pt_indices.size() == nb_pts) {
                auto cvx_type = M_params.domain.meshType == "GT_PK(3,1)"? 10 : 12;
                cells.push_back(reorderNodesForVTK(pt_indices, cvx_type == 12));
                cell_types.push_back(cvx_type);
            }
        }

        // Export to vtk (custom)
        std::string filename = "result_" + std::to_string(i) + ".vtk";

        std::ofstream out(filename);
        if (!out) {
            std::cerr << "Cannot open file " << filename << " for writing.\n";
            return;
        }

        out << "# vtk DataFile Version 3.0\n";
        out << "Custom GetFEM Export\n";
        out << "ASCII\n";
        out << "DATASET UNSTRUCTURED_GRID\n";

        out << "POINTS " << dof_coords.size() << " float\n";
        for (const auto& pt : dof_coords) {
            for (size_t d = 0; d < dim; ++d)
                out << pt[d] << " ";
            out << "\n";
        }

        out << "CELLS " << cells.size() << " ";
        size_type total_idx_count = 0;
        for (const auto &c : cells) total_idx_count += c.size() + 1;
        out << total_idx_count << "\n";

        for (const auto &c : cells) {
            out << c.size();
            for (auto idx : c) out << " " << idx;
            out << "\n";
        }

        out << "CELL_TYPES " << cell_types.size() << "\n";
        for (auto t : cell_types) out << t << "\n";

        out << "POINT_DATA " << U_vector.size() << "\n";
        out << "VECTORS U float\n";
        for (const auto& v : U_vector) {
            for (size_t d = 0; d < dim; ++d)
                out << v[d] << " ";
            out << "\n";
        }

        out << "FIELD StressField 1\n";
        out << "StressVoigt 6 " << stress_voigt.size() << " float\n";
        for (const auto &v : stress_voigt)
            out << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << " " << v[5] << "\n";

        std::cout << "done." << std::endl;

        M_model.set_real_variable("uL") = uL_backup;
        M_model.set_real_variable("uR") = uR_backup;

    }

} // namespace gf




        // INTERNAL EXPORT

        // getfem::vtk_export exp("result_" + std::to_string(i) + ".vtk");
        // exp.exporting(M_FEM.mf_u1());
        // exp.write_mesh();
        // exp.write_point_data(M_FEM.mf_u1(), M_model.real_variable("uL"), "uL");
        // exp.write_point_data(M_FEM.mf_u2(), M_model.real_variable("uR"), "uR");

        // gmm::clean(VL, 1E-20);
        // exp.exporting(M_FEM.mf_stress());
        // exp.write_point_data(M_FEM.mf_stress(), VL, "stressL");

        // gmm::clean(VR, 1E-20);
        // exp.exporting(M_FEM.mf_stress());
        // exp.write_point_data(M_FEM.mf_stress(), VR, "stressR");

        // UNCOMMENT THIS!!!
        // M_model.set_real_variable("uL") = uL_backup;
        // M_model.set_real_variable("uR") = uR_backup;

            // Alternatively, export in Voigt notation
                // size_type voigt_dim = (dim == 2) ? 3 : 6;
                // M_FEM.mf_stress().set_qdim(voigt_dim);

                // plain_vector VL(M_FEM.mf_stress().nb_dof());
                // plain_vector VR(M_FEM.mf_stress().nb_dof());

                // getfem::ga_interpolation_Lagrange_fem(M_model, "stressL_voigt", M_FEM.mf_stress(), VL);
                // getfem::ga_interpolation_Lagrange_fem(M_model, "stressR_voigt", M_FEM.mf_stress(), VR);
   