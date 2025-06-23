#ifndef _CONTACT_ENFORCEMENT_STRATEGY_HPP_
#define _CONTACT_ENFORCEMENT_STRATEGY_HPP_

#include "Core.hpp"
#include "Params.hpp"

namespace gf {

    /**
     * @brief Base class for contact enforcement strategies.
     * This class defines the interface for enforcing contact conditions in a finite element context.
     */
    class ContactEnforcementStrategy {
    public:
        virtual ~ContactEnforcementStrategy() = default;

        /**
         * @brief Enforce contact conditions on the given mesh_fem.
         * @param context The context containing the necessary data for enforcement.
         */
        virtual void enforce(getfem::model& md, const getfem::mesh_im& im) const = 0;
    };


    /**
     * @brief Nitsche contact enforcement strategy.
     * This class implements the Nitsche method for enforcing contact conditions.
     */
    class NitscheContactEnforcement : public ContactEnforcementStrategy {
        
        scalar_type M_theta; ///< Nitsche parameter theta
        scalar_type M_gammaN; ///< Nitsche parameter gamma0

    public:

        /**
         * @brief Constructor that initializes the Nitsche parameters.
         * @param theta Nitsche parameter theta
         * @param gammaN Nitsche parameter gamma0
         */
        NitscheContactEnforcement(scalar_type theta, scalar_type gammaN)
        : M_theta(theta), M_gammaN(gammaN) {}

        void enforce(getfem::model& md, const getfem::mesh_im& im) const override;
    };


    /**
     * @brief Penalty contact enforcement strategy.
     * This class implements the penalty method for enforcing contact conditions.
     * The penalty method applies a penalty term to the contact conditions.
     */
    class PenaltyContactEnforcement : public ContactEnforcementStrategy {

        scalar_type M_gammaP; ///< Penalty parameter gammaP

    public:

        /**
         * @brief Constructor for the PenaltyContactEnforcement class.
         * @param gammaP Penalty parameter gammaP
         */
        PenaltyContactEnforcement(scalar_type gammaP)
        : M_gammaP(gammaP) {}
    

        void enforce(getfem::model& md, const getfem::mesh_im& im) const override;
    };


    /**
     * @brief Augmented Lagrangian contact enforcement strategy.
     * This class implements the Augmented Lagrangian method for enforcing contact conditions.
     * It uses a parameter gammaL to control the enforcement strength.
     */
    class AugmentedLagrangianContactEnforcement : public ContactEnforcementStrategy {
        getfem::mesh_fem& M_mfLMn; ///< Mesh finite element for normal Lagrange multipliers
        getfem::mesh_fem& M_mfLMt; ///< Mesh finite element for tangential Lagrange multipliers
        scalar_type M_gammaL; ///< Augmented Lagrangian parameter gammaL

    public:
        /**
         * @brief Constructor for the AugmentedLagrangianContactEnforcement class.
         * @param gammaL Augmented Lagrangian parameter gammaL
         */
        AugmentedLagrangianContactEnforcement(
            scalar_type gammaL,
            getfem::mesh_fem& mfLMn,
            getfem::mesh_fem& mfLMt
        ): M_gammaL(gammaL), M_mfLMn(mfLMn), M_mfLMt(mfLMt){}

        void enforce(getfem::model& md, const getfem::mesh_im& im) const override;

    };


} // namespace gf



#endif // _CONTACT_ENFORCEMENT_STRATEGY_HPP_