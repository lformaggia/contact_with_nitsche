// #ifndef _KKT_BRICK_HPP_
// #define _KKT_BRICK_HPP_

// #include "Core.hpp"

// namespace gf {

//     size_type
//     add_Nitsche_KKT_condition_brick(
//         getfem::model &md,
//         const getfem::mesh_im &mim,
//         const std::string &varname_u1, const std::string &varname_u2,
//         const std::string &dataname_gamma0_coeff, const std::string &dataname_theta_coeff,
//         // size_type region);
//     struct Nitsche_KKT_condition_brick : public getfem::virtual_brick {

//         void asm_real_tangent_terms(const getfem::model& md, size_type ib,
//                                             const varnamelist& varl,
//                                             const varnamelist& datal,
//                                             const getfem::model::mimlist& mims,
//                                             getfem::model::real_matlist &matl,
//                                             getfem::model::real_veclist &vecl,
//                                             getfem::model::real_veclist &vecl_sym,
//                                             size_type region, build_version nl) const override;

//         Nitsche_KKT_condition_brick(scalar_type theta)
//         {
//             set_flags("Nitsche KKT brick", false, /* linear*/
//                                            theta == 1, /* symmetric*/
//                                            false, /* coercivity */
//                                            true, /* real version defined */
//                                            false /* no complex version */);
//         }
//     };

//     template<typename MAT, typename VECT>
//     void 
//     asm_Nitsche_KKT_nonlinear_matrix(MAT& M,
//                                      const getfem::mesh_im& mim,
//                                      const getfem::mesh_fem& mf,
//                                      const getfem::mesh_region &rg);


// } // namespace gf




// #endif // _KKT_BRICK_HPP