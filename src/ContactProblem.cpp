#include "ContactProblem.hpp"

namespace gf{

    bool contact_flag = false;

    ContactProblem::ContactProblem(const Mesh& m, const Params& p)
    : M_mesh(m),
    M_params(p),
    M_BC(m.get()),
    M_FEM(m.get()),
    M_integrationMethod(m.get()),
    M_integrationMethodSurface(m.get()),
    M_imData(M_integrationMethod)
    {
    }

    void
    ContactProblem::init()
    {

        M_BC.readBC(M_params.datafile);

        M_FEM.setMeshFem(M_params.datafile, M_mesh.get());
        

        getfem::pintegration_method ppi = getfem::int_method_descriptor(M_params.numerics.integration);
        size_type N = M_params.domain.dim;
        M_integrationMethod.set_integration_method(M_mesh.get().convex_index(), ppi);
        M_imData.set_tensor_size(bgeot::multi_index(N,N)); /** \todo */

        // getfem::pintegration_method ppiSurf = getfem::int_method_descriptor("IM_QUAD(5)");
        // M_integrationMethodSurface.set_integration_method(M_mesh.region(Fault).index(), ppiSurf);

    }


    void
    ContactProblem::assemble()
    {
        
        gmm::set_traces_level(1);

        dim_type dim = M_mesh.get().dim();
        size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();

        // Main unknown of the problem (displacement)
        M_model.add_fem_variable("uL", M_FEM.mf_u1());
        M_model.add_fem_variable("uR", M_FEM.mf_u2());

        // Add scalar data to the model
        /** \todo Change using the function element_size() from getfem */
        scalar_type h = M_params.domain.Lx / M_params.domain.Nx;  // typical mesh size

        M_model.add_initialized_scalar_data("lambda", M_params.physics.M_lambda);
        M_model.add_initialized_scalar_data("mu", M_params.physics.M_mu);
        M_model.add_initialized_scalar_data("gammaN", 10*M_params.physics.M_E0/h);
        M_model.add_initialized_scalar_data("theta", M_params.nitsche.theta);  // symmetric variant
        M_model.add_initialized_scalar_data("mu_fric", M_params.physics.M_mu_friction);


        // Add isotropic elasticity bricks
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uL", "lambda", "mu", RegionType::BulkLeft);
        getfem::add_isotropic_linearized_elasticity_brick(
            M_model, M_integrationMethod, "uR", "lambda", "mu", RegionType::BulkRight);
        std::cout << "Added elasticity bricks." << std::endl;
       

        // Define some useful macro for Nitsche contact integrals
        M_model.add_initialized_scalar_data("eps", 1.e-20);
        M_model.add_macro("n", "Normal"); // use normal on contact face
        /** \todo: add macros for t1, t2 */
        // Choose auxiliary vector `a` robustly
M_model.add_macro("eps_dir", "1e-4"); // to avoid numerical degeneracies
M_model.add_macro("a", "[0, 1, 0]");

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

        M_model.add_macro("traction_u", "((lambda*Trace(Grad_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_uL)) * n)");
        M_model.add_macro("traction_v", "((lambda*Trace(Grad_Test_uL)*Id(qdim(uL)) + 2*mu*Sym(Grad_Test_uL)) * n)");
        M_model.add_macro("sig_u_nL", "traction_u . n");
        M_model.add_macro("sig_u_t1", "traction_u . t1");
        M_model.add_macro("sig_u_t2", "traction_u . t2");
        M_model.add_macro("sig_v_nL", "traction_v . n");
        M_model.add_macro("sig_v_t1", "traction_v . t1");
        M_model.add_macro("sig_v_t2", "traction_v . t2");

        // Normal gap and stress
        M_model.add_macro("Pn_u", "(gammaN * un_jump - sig_u_nL)");
        M_model.add_macro("Pt1_u", "(gammaN * ut1_jump - sig_u_t1)");
        M_model.add_macro("Pt2_u", "(gammaN * ut2_jump - sig_u_t2)");
        M_model.add_macro("Pn_v_theta", "(gammaN * vn_jump - theta*sig_v_nL)");
        M_model.add_macro("Pt1_v_theta", "(gammaN * vt1_jump - theta*sig_v_t1)");
        M_model.add_macro("Pt2_v_theta", "(gammaN * vt2_jump - theta*sig_v_t2)");
        
        // Friction threshold
        M_model.add_macro("Sh", "mu_fric * pos_part(Pn_u)");
        M_model.add_macro("norm_Pt", "sqrt(Pt1_u*Pt1_u + Pt2_u*Pt2_u)");
        M_model.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
        M_model.add_macro("proj_Pt2_u", "Pt2_u * min(1, Sh / (norm_Pt + eps))");

        std::cout << "Added macros." << std::endl;


        // Add linear stress brick
        getfem::add_linear_term(
            M_model,
            M_integrationMethod,
            "- theta/gammaN * sig_u_nL * sig_v_nL", /** expression */
            Fault, /** region */
            false, /** symmetric */
            false, /** coercive */
            "linear_stress",
            false /** check */
        );


        // Add KKT condition brick
        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "1/gammaN * pos_part(Pn_u) * Pn_v_theta",
            Fault,
            false,
            false,
            "KKTbrick"
        );

        // Add Coulomb condition brick
        getfem::add_nonlinear_term(
            M_model,
            M_integrationMethod,
            "(1/gammaN) * (proj_Pt1_u * Pt1_v_theta + proj_Pt2_u * Pt2_v_theta)",
            Fault,
            false,
            false,
            "CoulombBrick"
        );


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
        plain_vector G(M_FEM.mf_rhs().nb_dof()*dim);

        for (size_type i = 0; i < nb_dof_rhs; ++i)
            gmm::copy(M_params.physics.M_gravity, gmm::sub_vector(G, gmm::sub_interval(i*dim, dim)));

        M_model.add_initialized_fem_data("VolumicData", M_FEM.mf_rhs(), G);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uL", "VolumicData", BulkLeft);
        getfem::add_source_term_brick(M_model, M_integrationMethod, "uR", "VolumicData", BulkRight);

        // Neumann conditions
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

        // Dirichlet conditions
        const auto & DirichletBCs = M_BC.Dirichlet();
        plain_vector D(nb_dof_rhs*dim);
        std::cout << "DirichletBCs.size(): " << DirichletBCs.size() << std::endl;

        for (const auto& bc: DirichletBCs){

            std::cout << "Assembling Dirichlet on boundary " << bc->ID() << std::endl;

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

        // Mixed conditions (normal Dirichlet)
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


        // Solve the problems
        gmm::iteration iter(M_params.it.atol, 1, M_params.it.maxIt);
        getfem::standard_solve(M_model,iter);

        // Export results
        getfem::vtk_export exp("result.vtk");
        exp.exporting(M_FEM.mf_u1());
        exp.write_mesh();
        exp.write_point_data(M_FEM.mf_u1(), M_model.real_variable("uL"), "uL");
        exp.write_point_data(M_FEM.mf_u2(), M_model.real_variable("uR"), "uR");

    }


} // namespace gf

        // const size_type Nb_t = 19;
        // std::vector<scalar_type> T(19);
        // T[0] = 0; T[1] = 0.9032; T[2] = 1; T[3] = 1.1; T[4] = 1.3;
        // T[5] = 1.5; T[6] = 1.7; T[7] = 1.74; T[8] = 1.7; T[9] = 1.5;
        // T[10] = 1.3; T[11] = 1.1; T[12] = 1; T[13] = 0.9032; T[14] = 0.7;
        // T[15] = 0.5; T[16] = 0.3; T[17] = 0.1; T[18] = 0;


        // for (size_type nb = 0; nb < Nb_t; ++nb) {
        //     cout << "=============iteration number : " << nb << "==========" << endl;
        
        //     scalar_type t = T[nb];
            
        //     // Defining the Neumann condition right hand side.
        //     base_small_vector v(dim);
        //     v[dim-1] = -2.0;
        //     gmm::scale(v,t);
            
        //     for (size_type i = 0; i < nb_dof_rhs; ++i)
        //     gmm::copy(v, gmm::sub_vector
        //         (F, gmm::sub_interval(i*dim, dim)));
            
        //     gmm::copy(F, M_model.set_real_variable("NeumannData"));
            
        //     // Generic solve.
        //     cout << "Number of variables : " 
        //     << M_model.nb_dof() << endl;
            
        //     getfem::newton_search_with_step_control ls;
        //     // getfem::simplest_newton_line_search ls;
        //     gmm::iteration iter(residual, 2, 40000);
        //     getfem::standard_solve(M_model, iter,
        //         getfem::rselect_linear_solver(M_model, "superlu"), ls);
        
            
            
        //     // Get the solution and save it
        //     gmm::copy(M_model.real_variable("u"), U);
            
            
        //     std::stringstream fname; fname << datafilename << "_" << nb << ".vtk";

        //     if (do_export) {
        //         getfem::vtk_export exp(fname.str());
        //         exp.write_point_data(M_FEM.mf_U1(), U, "displacement");
        //     }
            
        // }


        

        // if (contact_flag){
        //     // Add the Nitsche contact brick with raytracing
        //     size_type brick_id = getfem::add_Nitsche_large_sliding_contact_brick_raytracing(
        //         M_model,
        //         true,                   // unbiased
        //         "gamma_contact",        // Nitsche parameter
        //         1e-6,                   // release distance
        //         "mu_friction",          // friction coefficient
        //         "alpha",                // alpha
        //         false,                  // symmetric version
        //         false                   // frame indifferent
        //     );

        //     // Add contact boundaries (using both master and slave sides for unbiased)
        //     getfem::add_contact_boundary_to_Nitsche_large_sliding_contact_brick(
        //         M_model, brick_id, M_integrationMethod,
        //         M_mesh.getRegions().at("Fault")->ID(),
        //         false, true, true, "uL"
        //     );
        //     getfem::add_contact_boundary_to_Nitsche_large_sliding_contact_brick(
        //         M_model, brick_id, M_integrationMethod,
        //         M_mesh.getRegions().at("Fault")->ID(),
        //         true, false, true, "uR"
        //     );
        // }    

        // // // Non-penetrability KKT condition brick
        // // M_model.add_initialized_scalar_data("theta", M_params.nitsche.theta);
        // // M_model.add_initialized_scalar_data("gamma0", M_params.nitsche.gamma0);

        // // add_Nitsche_KKT_brick("u", /* faultRegionID */);

        // // // Coulomb friction condition brick
        // // add_Nitsche_friction_brick("u", /* faultRegionID */);

    

        // // // Dirichlet condition 
        // // /** \todo */




        // M_model.add_fem_data("stressL", M_FEM.mf_stress1());
        // M_model.add_fem_data("stressR", M_FEM.mf_stress2());
        // M_model.add_fem_data("stressL0", M_FEM.mf_stress1());
        // M_model.add_fem_data("stressR0", M_FEM.mf_stress2());
        // /** \todo: add other fem data */





        // size_type nb_dof_rhs = M_FEM.mf_rhs().nb_dof();
        // size_type N = M_mesh.dim();

        // Displacement: 0 is for the previous solution
        /**
         * \todo: alternatively, try to build a unique mesh_fem allover the mesh
         *    model.add_filtered_fem_variable("uL", M_FEM.mf_u(), BulkLeft.getRegion());
         */



        // model.add_im_data("Previous_Ep", mim_data);

        /* choose the projection type */
        // getfem::pconstraints_projection
        //     proj = std::make_shared<getfem::VM_projection>(0);

        // std::vector<std::string> plastic_variables = {"u", "xi", "Previous_Ep"};
        // std::vector<std::string> plastic_data = {"lambda", "mu", "sigma_y"};
        

        // add_small_strain_elastoplasticity_brick
        //     (model, mim, "Prandtl Reuss", getfem::DISPLACEMENT_ONLY,
        //     plastic_variables, plastic_data);
        
        // plain_vector F(nb_dof_rhs * N);
        // model.add_initialized_fem_data("NeumannData", M_FEM.mf_rhs(), F);


        // getfem::add_source_term_brick
        //     (model, mim, "u", "NeumannData", NEUMANN_BOUNDARY_NUM);
        
        // model.add_initialized_fem_data("DirichletData", mf_rhs, F);
        // getfem::add_Dirichlet_condition_with_multipliers
        //     (model, mim, "u", mf_u, DIRICHLET_BOUNDARY_NUM, 
        //     "DirichletData");

        // const size_type Nb_t = 19;
        // std::vector<scalar_type> T(19);
        // T[0] = 0; T[1] = 0.9032; T[2] = 1; T[3] = 1.1; T[4] = 1.3;
        // T[5] = 1.5; T[6] = 1.7; T[7] = 1.74; T[8] = 1.7; T[9] = 1.5;
        // T[10] = 1.3; T[11] = 1.1; T[12] = 1; T[13] = 0.9032; T[14] = 0.7;
        // T[15] = 0.5; T[16] = 0.3; T[17] = 0.1; T[18] = 0;


        // getfem::mesh_fem mf_vm(mesh);
        // mf_vm.set_classical_discontinuous_finite_element(1);
        // getfem::base_vector VM(mf_vm.nb_dof());
        // getfem::base_vector plast(mf_vm.nb_dof());
        
        // for (size_type nb = 0; nb < Nb_t; ++nb) {
        //     cout << "=============iteration number : " << nb << "==========" << endl;
        
        //     scalar_type t = T[nb];
            
        //     // Defining the Neumann condition right hand side.
        //     base_small_vector v(N);
        //     v[N-1] = -PARAM.real_value("FORCE");
        //     gmm::scale(v,t); 
            
        //     for (size_type i = 0; i < nb_dof_rhs; ++i)
        //     gmm::copy(v, gmm::sub_vector
        //         (F, gmm::sub_interval(i*N, N)));
            
        //     gmm::copy(F, model.set_real_variable("NeumannData"));
            
        //     // Generic solve.
        //     cout << "Number of variables : " 
        //     << model.nb_dof() << endl;



    // size_type
    // ContactProblem::add_Nitsche_KKT_brick(const varnamelist& u, const varnamelist& stress, size_type region)
    // {
    //     getfem::pbrick pbr = std::make_shared<my_KKT_Nitsche_brick>();
    //     getfem::model::termlist tl;
    //     tl.push_back(getfem::model::term_description(u[0], u[1], false));

    //     return M_model.add_brick(pbr, u, stress, tl,
    //                              getfem::model::mimlist(1, &M_integrationMethod),
    //                              region);
    // }

    // size_type
    // ContactProblem::add_Nitsche_friction_brick(const varnamelist& u, const varnamelist& stress, size_type region)
    // {
    //     getfem::pbrick pbr = std::make_shared<my_friction_Nitsche_brick>();
    //     getfem::model::termlist tl;
    //     tl.push_back(getfem::model::term_description(u[0], u[1], false));

    //     return M_model.add_brick(pbr, u, stress, tl,
    //                              getfem::model::mimlist(1, &M_integrationMethod),
    //                              region);
    // }

