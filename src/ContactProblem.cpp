#include "ContactProblem.hpp"

namespace gf{

    bool contact_flag = false;

    ContactProblem::ContactProblem(const Mesh& m, const Params& p)
    : M_mesh(m),
    M_params(p),
    M_BC(m.get()),
    M_FEM(m.get()),
    M_integrationMethod(m.get())
    {
    }

    void
    ContactProblem::init()
    {

        std::cout << "Initializing the contact problem...";
        M_BC.readBC(M_params.datafile);
        
        M_FEM.setMeshFem(M_params.numerics);

        getfem::pintegration_method ppi = getfem::int_method_descriptor(M_params.numerics.integration);
        size_type N = M_params.domain.dim;
        M_integrationMethod.set_integration_method(M_mesh.get().convex_index(), ppi);
    
        // Select the method    
        if (M_params.contact.method == "nitsche")
            M_contactEnforcement = std::make_unique<NitscheContactEnforcement>(
                M_params.contact.theta, M_params.contact.gammaN);
        else if (M_params.contact.method == "penalty")
            M_contactEnforcement = std::make_unique<PenaltyContactEnforcement>(
                M_params.contact.gammaP);
        else if (M_params.contact.method == "augLM")
            M_contactEnforcement = std::make_unique<AugmentedLagrangianContactEnforcement>(
                M_params.contact.gammaL, M_FEM.mf_LM());
        else
            throw std::runtime_error("Unknown contact enforcement method: " + M_params.contact.method);
    
        std::cout << "done." << std::endl;
        
    }


    void
    ContactProblem::assemble()
    {
        gmm::set_traces_level(M_params.verbose ? 1 : 0);

        std::cout << "Preparing the assembly phase..." << std::endl;
        
        dim_type dim = M_mesh.get().dim();
        size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();

        // Main unknown of the problem (displacement)
        if (M_params.verbose) std::cout << "  Defining variables...";
        M_model.add_fem_variable("uL", M_FEM.mf_u());
        M_model.add_fem_variable("uR", M_FEM.mf_u());
        if (M_params.verbose) std::cout << "done.\n";

        // Add scalar data to the model
        if (M_params.verbose) std::cout << "  Initializing scalar data...";
        M_model.add_initialized_scalar_data("lambda", M_params.physics.M_lambda);
        M_model.add_initialized_scalar_data("mu", M_params.physics.M_mu);
        M_model.add_initialized_scalar_data("mu_fric", M_params.physics.M_mu_friction);
        if (M_params.verbose) std::cout << "done.\n";

        // Define some useful macro for Nitsche contact integrals
        if (M_params.verbose) std::cout << "  Defining macros...";
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

        // Enforcement of contact conditions
        M_contactEnforcement->enforce(M_model, M_integrationMethod, M_params.verbose);

        // Add isotropic elasticity bricks
        if (M_params.verbose) std::cout << "  Adding elasticity bricks...";
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uL", "lambda", "mu", RegionType::BulkLeft);
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uR", "lambda", "mu", RegionType::BulkRight);
        if (M_params.verbose) std::cout << "done.\n";

        // Stabilizers for displacement variables (negligible, eps = 1.e-20)
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
        if (M_params.verbose) std::cout << "  Adding volumic source term brick...";
        plain_vector G(M_FEM.mf_rhs().nb_dof()*dim);

        for (size_type i = 0; i < nb_dof_rhs; ++i)
            gmm::copy(M_params.physics.M_gravity, gmm::sub_vector(G, gmm::sub_interval(i*dim, dim)));

        M_model.add_initialized_fem_data("VolumicData", M_FEM.mf_rhs(), G);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uL", "VolumicData", BulkLeft);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uR", "VolumicData", BulkRight); 
        if (M_params.verbose) std::cout << "done.\n";

        // Neumann conditions
        if (M_params.verbose) std::cout << "  Adding Neumann condition bricks...";
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
        if (M_params.verbose) std::cout << "done.\n";

        // Dirichlet conditions
        if (M_params.verbose) std::cout << "  Adding Dirichlet condition bricks...";
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
                    (M_model, M_integrationMethod, "uL", M_FEM.mf_u(), bc->ID(), 
                     bc->name());
            else
                getfem::add_Dirichlet_condition_with_multipliers
                    (M_model, M_integrationMethod, "uR", M_FEM.mf_u(), bc->ID(), 
                     bc->name());
        }
        if (M_params.verbose) std::cout << "done.\n";

        // Mixed conditions (normal Dirichlet)
        if (M_params.verbose) std::cout << "  Adding normal Dirichlet condition bricks...";
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
        if (M_params.verbose) std::cout << "done.\n";

    }

    void
    ContactProblem::solve()
    {   
        plain_vector uL0 (M_FEM.mf_u().nb_basic_dof(),0);
        plain_vector uR0 (M_FEM.mf_u().nb_basic_dof(),0);

        if (M_params.init){
            // the input files are used to better initialize the solution at the first time-step
            importCSV(uL0, 0, "initSolutionLeft.csv");
            importCSV(uR0, 0, "initSolutionRight.csv");
            M_model.set_real_variable("uL") = uL0;
            M_model.set_real_variable("uR") = uR0;
        }

        if (M_params.verbose) std::cout << "Solving the problem..." << std::endl;

        const auto & NeumannBCs = M_BC.Neumann();
        
        dim_type dim = M_mesh.get().dim();
        size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();
        plain_vector F(nb_dof_rhs*dim);
        
        plain_vector VM(M_FEM.mf_stress().nb_dof());
        
        // Time loop
        size_type n_timesteps = static_cast<size_type>(M_params.time.tend - M_params.time.t0)/M_params.time.dt;
        
        scalar_type t {};

        for (size_type i{}; i < n_timesteps; ++i)
        {
            std::cout << "t = " << t << std::endl;
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
            int noise_level = 0;
            if (M_params.verbose) noise_level = 1;
            gmm::iteration iter(M_params.it.tol, noise_level, M_params.it.maxIt);

            getfem::standard_solve(M_model,iter);

            if (M_params.verbose) {
                if (iter.converged()) {
                    std::cout << "  Iteration converged after " << iter.get_iteration() << " iterations." << std::endl;
                } else {
                    std::cerr << " Warning: Iteration did not converge after " << iter.get_iteration() << " iterations." << std::endl;
                }
            }

            /* Output solution to csv file:
                - if -r: refSolution{Left,Right}.csv will be used as reference for computing the erorr
                - if -ir: initSolution{Left,Right}.csv will be used to initialize uL, uR
            */
            if (M_params.refined || M_params.initRef){
                size_type dofU = M_FEM.mf_u().nb_basic_dof();
                // Export the refined solution
                std::vector<scalar_type> uL(dofU);
                gmm::copy(M_model.real_variable("uL"), uL);
                std::vector<scalar_type> uR(dofU);
                gmm::copy(M_model.real_variable("uR"), uR);
                std::string outfileL = (M_params.refined) ? "refSolutionLeft.csv" : "initSolutionLeft.csv";
                std::string outfileR = (M_params.refined) ? "refSolutionRight.csv" : "initSolutionRight.csv";
                exportCSV(uL, i, outfileL);
                exportCSV(uR, i, outfileR);
            }
            else if (M_params.test) {
                computeError(i);
            }

            // Export results
            exportVtk(i);

            // Advance in time
            t += M_params.time.dt;
        }

    }


    void
    ContactProblem::computeError(size_type i) {
        // Create a finer mesh for the "exact" solution (better do this elsewhere)
        Params p(M_params);
        // RMK: p.domain.h 
        p.domain.h = 0.125;
        int ratio = static_cast<int>(p.domain.Ly/p.domain.Lx);
        p.domain.Nx = static_cast<int>(p.domain.Lx/p.domain.h);
        p.domain.Ny = static_cast<int>(p.domain.Ly/p.domain.h/ratio);
        p.domain.Nz = static_cast<int>(p.domain.Lz/p.domain.h/ratio);
        Mesh m(p);
        getfem::mesh_fem mfUrefined(m.get(),3);
        mfUrefined.set_finite_element(getfem::fem_descriptor(M_params.numerics.FEMTypeDisplacement));
        size_type dofUrefined = mfUrefined.nb_basic_dof();
        size_type dofU = M_FEM.mf_u().nb_basic_dof();
        std::vector<scalar_type> uL_refined(dofUrefined);
        std::vector<scalar_type> uR_refined(dofUrefined);

        importCSV(uL_refined, i, "refSolutionLeft.csv");
        importCSV(uR_refined, i, "refSolutionRight.csv");

        std::vector<scalar_type> uL(dofU);
        gmm::copy(M_model.real_variable("uL"), uL);
        std::vector<scalar_type> uR(dofU);
        gmm::copy(M_model.real_variable("uR"), uR);

        // Compute error w.r.t the refined solution
        scalar_type errL2 {};
        scalar_type errH1 {};

        plain_vector uL_refined_interp(dofU,0.);
        plain_vector uR_refined_interp(dofU,0.);

        getfem::interpolation(mfUrefined, M_FEM.mf_u(), uL_refined, uL_refined_interp);
        getfem::interpolation(mfUrefined, M_FEM.mf_u(), uR_refined, uL_refined_interp);
        errL2 += getfem::asm_L2_dist(M_integrationMethod, M_FEM.mf_u(), uL, M_FEM.mf_u(), uL_refined_interp, M_mesh.region(BulkLeft));
        errH1 += getfem::asm_H1_dist(M_integrationMethod, M_FEM.mf_u(), uL, M_FEM.mf_u(), uL_refined_interp, M_mesh.region(BulkLeft));
        errL2 += getfem::asm_L2_dist(M_integrationMethod, M_FEM.mf_u(), uR, M_FEM.mf_u(), uR_refined_interp, M_mesh.region(BulkRight));
        errH1 += getfem::asm_H1_dist(M_integrationMethod, M_FEM.mf_u(), uR, M_FEM.mf_u(), uR_refined_interp, M_mesh.region(BulkRight));


        
        // getfem::interpolation(M_FEM.mf_u(), mfUrefined, uL_refined_interp, uL_refined);
        // // getfem::interpolation(M_FEM.mf_u(), mfUrefined, uL_refined_interp, uR_refined);
        // errL2 += getfem::asm_L2_dist(M_integrationMethod, M_FEM.mf_u(), uL_refined_interp, M_FEM.mf_u(), uL, M_mesh.region(BulkLeft));
        // errH1 += getfem::asm_H1_dist(M_integrationMethod, M_FEM.mf_u(), uL_refined_interp, M_FEM.mf_u(), uL, M_mesh.region(BulkLeft));
        // errL2 += getfem::asm_L2_dist(M_integrationMethod, M_FEM.mf_u(), uR_refined_interp, M_FEM.mf_u(), uR, M_mesh.region(BulkRight));
        // errH1 += getfem::asm_H1_dist(M_integrationMethod, M_FEM.mf_u(), uR_refined_interp, M_FEM.mf_u(), uR, M_mesh.region(BulkRight));

        std::cout << "  ErrL2: " << errL2 << std::endl;
        std::cout << "  ErrH1: " << errH1 << std::endl;
    }

    
    void ContactProblem::exportVtk(size_type i) {

        dim_type dim = M_mesh.get().dim();

        // Zero out the noised components
        plain_vector uL_backup = M_model.real_variable("uL");
        plain_vector uR_backup = M_model.real_variable("uR");

        plain_vector &uL = M_model.set_real_variable("uL");
        plain_vector &uR = M_model.set_real_variable("uR");

        auto dofs_uL = M_FEM.mf_u().dof_on_region(M_mesh.region(BulkLeft));
        auto dofs_uR = M_FEM.mf_u().dof_on_region(M_mesh.region(BulkRight));
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
        if (M_params.verbose) std::cout << "Computing stress...";

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

        if (M_params.verbose) std::cout << "done." << std::endl;

        // Custom export to VTK
        if (M_params.verbose) std::cout << "Exporting results to result_" << i << ".vtk...";

        scalar_type offset = 1e-5; // offset for visualization

        // const plain_vector& uL = M_model.real_variable("uL");
        // const plain_vector& uR = M_model.real_variable("uR");

        size_type nb_dof_L = dofs_uL.card();
        size_type nb_dof_R = dofs_uR.card();
        size_type nb_dof_F = dofs_uF.card();
        size_type nb_dof = M_FEM.mf_u().nb_dof();
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
                base_node pt = M_FEM.mf_u().point_of_basic_dof(dof);
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
                base_node pt = M_FEM.mf_u().point_of_basic_dof(dof);
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
        std::string filename = "output/result_" + std::to_string(i) + ".vtk";

        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Cannot open " + filename + " for writing.\n");
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

        if (M_params.verbose) std::cout << "done." << std::endl;

        M_model.set_real_variable("uL") = uL_backup;
        M_model.set_real_variable("uR") = uR_backup;

    }


    void
    ContactProblem::exportCSV(const plain_vector& U, size_type t, const std::string& filename)
    {
        std::ofstream file;
        if (t == 0) file.open(filename); // start writing
        else file.open(filename, std::ios::app); // open in append mode for next timesteps
        if (!file.is_open()) {
            throw std::runtime_error("Error: cannot open file " + filename + " for writing.\n");
        }

        const std::size_t n_nodes = U.size()/3;
        file << "t = " << t << "\n";

        size_type dof_per_node = 3;
        for (size_type n = 0; n < n_nodes; ++n) {
            for (size_type d = 0; d < 3; ++d) {
                file << std::setprecision(20) << U[n * dof_per_node + d];
                if (d < dof_per_node - 1)
                    file << ",";
            }
            file << "\n";
        }

        file << "\n";
        file.close();

    }


    void
    ContactProblem::importCSV(plain_vector& U, size_type i, const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening file " + filename + "\n");
        }

        std::string line;
        bool found_time_block = false;
        std::ostringstream time_marker;
        time_marker << "t = " << i;

        // Search for the line "t = i"
        while (std::getline(file, line)) {
            if (line.find(time_marker.str()) != std::string::npos) {
                found_time_block = true;
                break;
            }
        }

        if (!found_time_block) {
            std::cerr << "Time block t = " << i << " not found in file.\n";
            return;
        }

        U.clear(); // Clear U before filling

        // Read lines until an empty line or EOF
        while (std::getline(file, line)) {
            if (line.empty()) break; // Stop at empty line

            std::istringstream ss(line);
            std::string value_str;
            while (std::getline(ss, value_str, ',')) {
                std::istringstream vs(value_str);
                scalar_type value;
                vs >> value;
                if (!vs.fail()) {
                    U.push_back(value);
                } else {
                    // std::cerr << "Warning: couldn't parse value '" << value_str << "'\n";
                }
            }
        }
    }


} // namespace gf
