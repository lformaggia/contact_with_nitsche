// #include "KKTbrick.hpp"

// namespace gf {

// //     template<typename MAT, typename VECT>
// //     void 
// //     asm_Nitsche_KKT_nonlinear_matrix(
// //         MAT& M,
// //         const getfem::mesh_im& mim,
// //         const getfem::mesh_fem& mf,
// //         const getfem::mesh_region &rg)
// //     {
// //         getfem::ga_workspace workspace;
// //         gmm::sub_interval Iu(0, mf.nb_dof());
// //         base_vector U(mf.nb_dof());
// //         workspace.add_fem_variable("u", mf, Iu, U);
// //         workspace.add_expression(" /* TODO */", mim, rg);
// //         workspace.set_assembled_matrix(M);
// //         workspace.assembly(2);
// //     }
    
//     void
//     Nitsche_KKT_condition_brick::asm_real_tangent_terms(
//                                     const getfem::model& md, size_type ib,
//                                     const varnamelist& varl,
//                                     const varnamelist& datal,
//                                     const getfem::model::mimlist& mims,
//                                     getfem::model::real_matlist &matl,
//                                     getfem::model::real_veclist &vecl,
//                                     getfem::model::real_veclist &vecl_sym,
//                                     size_type rg, build_version nl) const {

//         GMM_ASSERT1(matl.size() == 1,
//                 "Nitsche KKT condition brick has one and only one term");
//         GMM_ASSERT1(mims.size() == 1,
//                 "Nitsche KKT condition brick needs one and only one mesh_im");
//         GMM_ASSERT1(varl.size() == 2,
//                 "Wrong number of variables for Nitsche KKT condition brick");
        
//         const getfem::mesh_fem &mf_u1 = md.mesh_fem_of_variable(varl[0]);
//         const getfem::mesh_fem &mf_u2 = md.mesh_fem_of_variable(varl[1]);
//         const getfem::mesh_im &mim = *mims[0];

//         gmm::clear(matl[0]);

//         getfem::mesh_region region = /*...*/;

//         asm_Nitsche_KKT_nonlinear_matrix(matl[0], mim, mf_u, region);
            
//     }

//     size_type
//     add_Nitsche_KKT_condition_brick(
//         getfem::model &md,
//         const getfem::mesh_im &mim,
//         const std::string &varname_u1, const std::string &varname_u2,
//         const std::string &dataname_gamma0_coeff, const std::string &dataname_theta_coeff,
//         size_type region)
//     {
//         scalar_type theta = md.real_variable(dataname_theta_coeff)[0];
//         getfem::pbrick pbr = std::make_shared<Nitsche_KKT_condition_brick>(theta);
//         getfem::model::termlist tl;
//         tl.push_back(getfem::model::term_description(varname_u1, varname_u2, theta==1));
//         return md.add_brick(pbr, varnamelist{varname_u1, varname_u2},
//                             varnamelist{dataname_gamma0_coeff,dataname_theta_coeff}, tl,
//                             getfem::model::mimlist(1, &mim), region);
//     }

// //     // std::string test_varname
// //     //   = "Test_" + sup_previous_and_dot_to_varname(varname);

// //     // std::string expr1 = "((("+dataexpr1+")*(Div_"+varname+"-Div_"+dataname3
// //     //   +"))*Id(meshdim)+(2*("+dataexpr2+"))*(Sym(Grad_"+varname
// //     //   +")-Sym(Grad_"+dataname3+"))):Grad_" +test_varname;
// //     // std::string expr2 = "(Div_"+varname+"*(("+dataexpr1+")*Id(meshdim))"
// //     //   +"+(2*("+dataexpr2+"))*Sym(Grad_"+varname+")):Grad_"+test_varname;

// //     // bool is_lin;
// //     // varnamelist vl, dl;
// //     // { // reenables disabled variables
// //     //   getfem::ga_workspace workspace(md, getfem::ga_workspace::inherit::ALL);
// //     //   workspace.add_expression(expr2, mim, region);
// //     //   varnamelist vl_test1, vl_test2;
// //     //   is_lin = workspace.used_variables(vl, vl_test1, vl_test2, dl, 2);
// //     // }
// //     //   return add_nonlinear_generic_assembly_brick
// //     //     (md, mim, dataname3.size() ? expr1 : expr2, region, false, false,
// //     //      "Linearized isotropic elasticity (with nonlinear dependance)");
// //     }


// } // namespace gf


// // struct Nitsche_large_sliding_contact_brick_raytracing
// //     : public virtual_brick {

// //     struct contact_boundary {
// //       size_type region;
// //       std::string varname_u;
// //       std::string sigma_u;
// //       std::string varname_w; // not needed, velocity
// //       bool is_master;
// //       bool is_slave;
// //       bool is_unbiased;
// //       const mesh_im *mim;
      
// //       std::string expr;
// //     };

// //     std::vector<contact_boundary> boundaries;
// //     std::string transformation_name;
// //     std::string u_group;
// //     std::string w_group; // velocity, not needed
// //     std::string friction_coeff;
// //     std::string alpha;
// //     std::string Nitsche_param; // gamma
// //     model::varnamelist vl, dl;
// //     model::mimlist ml;

// //     bool  sym_version, frame_indifferent, unbiased;

// //     void add_contact_boundary(model &md, const mesh_im &mim, size_type region,
// //                               bool is_master, bool is_slave, bool is_unbiased,
// //                               const std::string &u,
// //                               const std::string &sigma_u,
// //                               const std::string &w = "") {
// //       std::string test_u = "Test_" + sup_previous_and_dot_to_varname(u);
// //       std::string test_u_group = "Test_" + sup_previous_and_dot_to_varname(u_group);
// //       GMM_ASSERT1(is_slave || is_master, "The contact boundary should be "
// //                   "either master, slave or both");
// //       if (is_unbiased) {
// //                         GMM_ASSERT1((is_slave && is_master), "The contact boundary should be "
// //                           "both master and slave for the unbiased version");
// //                         is_slave=true; is_master=true;
// //                        }
// //       const mesh_fem *mf = md.pmesh_fem_of_variable(u);
// //       GMM_ASSERT1(mf, "The displacement variable should be a f.e.m. one");
// //       GMM_ASSERT1(&(mf->linked_mesh()) == &(mim.linked_mesh()),
// //                   "The displacement variable and the integration method "
// //                   "should share the same mesh");

// //       if (w.size()) {
// //         const mesh_fem *mf2 =  md.pmesh_fem_of_variable(w);
// //         GMM_ASSERT1(!mf2 || &(mf2->linked_mesh()) == &(mf->linked_mesh()),
// //                     "The data for the sliding velocity should be defined on "
// //                     " the same mesh as the displacement variable");
// //       }

// //       for (size_type i = 0; i < boundaries.size(); ++i) {
// //         const contact_boundary &cb = boundaries[i];
// //         if (&(md.mesh_fem_of_variable(cb.varname_u).linked_mesh())
// //             == &(mf->linked_mesh()) && cb.region == region)
// //           GMM_ASSERT1(false, "This contact boundary has already been added");
// //       }
// //       if (is_master)
// //         add_master_contact_boundary_to_raytracing_transformation
// //           (md, transformation_name, mf->linked_mesh(), u_group, region);
// //       else
// //          add_slave_contact_boundary_to_raytracing_transformation
// //           (md, transformation_name, mf->linked_mesh(), u_group, region);
      
// //       boundaries.push_back(contact_boundary());
// //       contact_boundary &cb = boundaries.back();
// //       cb.region = region;
// //       cb.varname_u = u;
// //       if (is_slave) cb.sigma_u = sigma_u;
// //       cb.varname_w = w;
// //       cb.is_master = is_master;
// //       cb.is_slave = is_slave;
// //       cb.is_unbiased= is_unbiased;
// //       cb.mim = &mim;
// //       if (is_slave) {
// //         std::string n, n0, Vs, g, Y, gamma;

// //         gamma ="("+Nitsche_param+"/element_size)";
// //         n = "Transformed_unit_vector(Grad_"+u+", Normal)";
// //         n0 = "Transformed_unit_vector(Grad_"+w+", Normal)";

// //         // For deformable bodies:
// //         // Coulomb_friction_coupled_projection(sigma(u),
// //         //    Transformed_unit_vector(Grad_u, Normal),
// //         //    (u-Interpolate(ug,trans)-(w-Interpolate(wg,trans)))*alpha,
// //         //    (Interpolate(X,trans)+Interpolate(ug,trans)-X-u).
// //         //      Transformed_unit_vector(Grad_u, Normal), f, r)
// //         Y = "Interpolate(X,"+transformation_name+")";
// //         g = "("+Y+"+Interpolate("+u_group+","+transformation_name+")-X-"+u+")."+n;
// //         Vs = "("+u+"-Interpolate("+u_group+","+transformation_name+")";
// //         if (w.size()) {
// //           Vs += "-"+w+"+Interpolate("+w_group+","+transformation_name+")";
// //           if (frame_indifferent)
// //             Vs += "+("+g+")*("+n+"-"+n0+")";
// //         }
// //         Vs += ")*"+alpha;

// //         std::string coupled_projection_def =
// //           "Coulomb_friction_coupled_projection("
// //           + sigma_u +","+n+","+Vs+","+g+","+friction_coeff+","
// //           + gamma +")";

// //         // For regid obstacle:
// //         // Coulomb_friction_coupled_projection(sigma(u),
// //         //   Transformed_unit_vector(Grad_u, Normal), (u-w)*alpha,
// //         //   (Interpolate(X,trans)-X-u).Transformed_unit_vector(Grad_u,
// //         //                                                      Normal), f, r)
// //         g = "(Interpolate(X,"+transformation_name+")-X-"+u+")."+n;
// //         if (frame_indifferent && w.size())
// //           Vs = "("+u+"-"+w+"+"+g+"*("+n+"-"+n0+"))*"+alpha;
// //         else
// //           Vs = "("+u+(w.size() ? ("-"+w):"")+")*"+alpha;

// //         std::string coupled_projection_rig =
// //           "Coulomb_friction_coupled_projection("
// //           + sigma_u +","+n+","+Vs+","+g+","+ friction_coeff+","
// //           + gamma +")";

// //         cb.expr =
// //           // 0.5* for non-biaised version
// //           (is_unbiased ? "-0.5*" : "-")
// //           // -coupled_projection_def.Test_u 
// //           + ("Interpolate_filter("+transformation_name+","
// //                             +coupled_projection_def+"."+test_u+",1) ") 
// //           +(is_unbiased ? "":"-Interpolate_filter("+transformation_name+","
// //                             +coupled_projection_rig+"."+test_u+",2) ")
// //           // Interpolate_filter(trans,
// //           //                   lambda.Interpolate(Test_ug, contact_trans), 1)
// //           // or
// //           // Interpolate_filter(trans,
// //           //     coupled_projection_def.Interpolate(Test_ug, contact_trans), 1)
// //           + (is_unbiased ? "+ 0.5*" : "+ ")
// //           +"Interpolate_filter("+transformation_name+","
// //                             +coupled_projection_def +".Interpolate("+test_u_group+"," + transformation_name+"), 1)"
// //           +(is_unbiased ? "":"+ Interpolate_filter("+transformation_name+","
// //                             +coupled_projection_rig +".Interpolate("+test_u_group+"," + transformation_name+"), 2)");                  
// //         }
// //       }

// //     virtual void asm_real_tangent_terms(const model &md, size_type ,
// //                                         const model::varnamelist &,
// //                                         const model::varnamelist &,
// //                                         const model::mimlist &,
// //                                         model::real_matlist &,
// //                                         model::real_veclist &,
// //                                         model::real_veclist &,
// //                                         size_type,
// //                                         build_version) const {
// //       // GMM_ASSERT1(mims.size() == 1,
// //       //            "Generic linear assembly brick needs one and only one "
// //       //            "mesh_im"); // to be verified ...

// //       for (const contact_boundary &cb : boundaries) {
// //         if (cb.is_slave)
// //           md.add_generic_expression(cb.expr, *(cb.mim), cb.region);
// //       }
// //     }


// //     Nitsche_large_sliding_contact_brick_raytracing
// //     ( bool unbia,
// //       const std::string &Nitsche_parameter,
// //      const std::string &f_coeff,const std::string &ug,
// //      const std::string &wg, const std::string &tr,
// //      const std::string &alpha_ = "1", bool sym_v = false,
// //      bool frame_indiff = false) {
// //       transformation_name = tr;
// //       u_group = ug; w_group = wg;
// //       friction_coeff = f_coeff;
// //       alpha = alpha_;
// //       Nitsche_param = Nitsche_parameter;
// //       sym_version = sym_v;
// //       frame_indifferent = frame_indiff;
// //       unbiased = unbia;

// //       set_flags("Integral large sliding contact bick raytracing",
// //                 false /* is linear*/,
// //                 false /* is symmetric */, false /* is coercive */,
// //                 true /* is real */, false /* is complex */);
// //     }

// //   };
