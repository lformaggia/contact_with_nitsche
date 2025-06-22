#include "ContactEnforcementStrategy.hpp"


namespace gf {

    void
    NitscheContactEnforcement::enforce(getfem::model& md, const getfem::mesh_im& im) const {
        
        md.add_initialized_scalar_data("gammaN", M_gamma0);
        md.add_initialized_scalar_data("theta", M_theta);  // symmetric variant
        
        // Normal gap and stress
        md.add_macro("Pn_u", "(gammaN * un_jump - sig_u_nL)");
        md.add_macro("Pt1_u", "(gammaN * ut1_jump - sig_u_t1)");
        md.add_macro("Pt2_u", "(gammaN * ut2_jump - sig_u_t2)");
        md.add_macro("Pn_v_theta", "(gammaN * vn_jump - theta*sig_v_nL)");
        md.add_macro("Pt1_v_theta", "(gammaN * vt1_jump - theta*sig_v_t1)");
        md.add_macro("Pt2_v_theta", "(gammaN * vt2_jump - theta*sig_v_t2)");

        // Friction threshold
        md.add_macro("Sh", "mu_fric * pos_part(Pn_u)");
        md.add_macro("norm_Pt", "sqrt(Pt1_u*Pt1_u + Pt2_u*Pt2_u)");
        md.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
        md.add_macro("proj_Pt2_u", "Pt2_u * min(1, Sh / (norm_Pt + eps))");

        std::cout << "done.\n";
        
        
        // Add linear stress brick
        std::cout << "  Adding linear stress brick...";
        getfem::add_linear_term(
            md,
            im,
            "- theta/gammaN * sig_u_nL * sig_v_nL", /** expression */
            Fault, /** region */
            false, /** symmetric */
            false, /** coercive */
            "linear_stress",
            false /** check */
        );
        std::cout << "done.\n";


        // Add KKT condition brick
        std::cout << "  Adding KKT condition brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "1/gammaN * pos_part(Pn_u) * Pn_v_theta",
            Fault,
            false,
            false,
            "KKTbrick"
        );
        std::cout << "done.\n";


        // Add Coulomb condition brick
        std::cout << "  Adding Coulomb friction brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "(1/gammaN) * (proj_Pt1_u * Pt1_v_theta + proj_Pt2_u * Pt2_v_theta)",
            Fault,
            false,
            false,
            "CoulombBrick"
        );
        std::cout << "done.\n";

    }



    void
    PenaltyContactEnforcement::enforce(getfem::model& md, const getfem::mesh_im& im) const {
        
        md.add_initialized_scalar_data("epsilon", M_epsilon);
        md.add_macro("Sh", "mu_fric * pos_part(un_jump)");
        md.add_macro("norm_ut_jump", "sqrt(ut1_jump*ut1_jump + ut2_jump*ut2_jump)");
        md.add_macro("proj_ut1_jump", "ut1_jump * min(1, Sh / (norm_ut_jump + eps))");
        md.add_macro("proj_ut2_jump", "ut2_jump* min(1, Sh / (norm_ut_jump + eps))");

        std::cout << "done.\n";
        
        // Add KKT condition brick
        std::cout << "  Adding KKT condition brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "epsilon * pos_part(un_jump) * vn_jump",
            Fault,
            false,
            false,
            "KKTbrick"
        );
        std::cout << "done.\n";


        // Add Coulomb condition brick
        std::cout << "  Adding Coulomb friction brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "epsilon*(proj_ut1_jump * vt1_jump + proj_ut2_jump * vt2_jump)",
            Fault,
            false,
            false,
            "CoulombBrick"
        );
        std::cout << "done.\n";

    }


    void
    AugmentedLagrangianContactEnforcement::enforce(getfem::model& md, const getfem::mesh_im& im) const {

        M_mfLM.set_qdim(1);
        md.add_fem_variable("lambdan", M_mfLM);
        M_mfLM.set_qdim(3);
        md.add_fem_variable("lambdat", M_mfLM);
        M_mfLM.set_qdim(1);
        
        md.add_initialized_scalar_data("gammaL", M_gammaL);
        
        // Normal gap and stress
        md.add_macro("lambdat1", "lambdat . t1");
        md.add_macro("lambdat2", "lambdat . t2");
        md.add_macro("mut1", "Test_lambdat . t1");
        md.add_macro("mut2", "Test_lambdat . t2");

        md.add_macro("Pn_u", "(gammaL * un_jump - lambdan)");
        md.add_macro("Pt1_u", "(gammaL * ut1_jump - lambdat1)");
        md.add_macro("Pt2_u", "(gammaL * ut2_jump - lambdat2)");

        // Friction threshold
        md.add_macro("Sh", "mu_fric * pos_part(Pn_u)");
        md.add_macro("norm_Pt", "sqrt(Pt1_u*Pt1_u + Pt2_u*Pt2_u)");
        md.add_macro("proj_Pt1_u", "Pt1_u * min(1, Sh / (norm_Pt + eps))");
        md.add_macro("proj_Pt2_u", "Pt2_u * min(1, Sh / (norm_Pt + eps))");

        std::cout << "done.\n";
        

        // Add KKT condition brick
        std::cout << "  Adding KKT condition brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "gammaL * pos_part(Pn_u) * vn_jump",
            Fault,
            false,
            false,
            "KKTbrick"
        );
        std::cout << "done.\n";

        M_mfLM.set_qdim(3);
        // Add Coulomb condition brick
        std::cout << "  Adding Coulomb friction brick...";
        getfem::add_nonlinear_term(
            md,
            im,
            "gammaL * (proj_Pt1_u * vt1_jump + proj_Pt2_u * vt2_jump)",
            Fault,
            false,
            false,
            "CoulombBrick"
        );
        std::cout << "done.\n";
        
        
        M_mfLM.set_qdim(1);
        // Add LM term for normal gap
        std::cout << "  Adding Lagrange multiplier term for normal gap...";
        getfem::add_nonlinear_term(
            md,
            im,
            " - 1/gammaL * pos_part(Pn_u) * Test_lambdan",
            Fault,
            false,
            false,
            "LM_NormalGapBrick"
        );
        std::cout << "done.\n";

        
        M_mfLM.set_qdim(1);
        // Add LM term for tangential gap
        std::cout << "  Adding Lagrange multiplier term for tangential gap...";
        getfem::add_nonlinear_term(
            md,
            im,
            " (lambdat1 + proj_Pt1) * mut1 + (lambdat2 + proj_Pt2) * mut2",
            Fault,
            false,
            false,
            "LM_TangentialGapBrick"
        );
        std::cout << "done.\n";
        
    }


} // namespace gf