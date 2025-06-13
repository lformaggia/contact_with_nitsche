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
        scalar_type M_gamma0; ///< Nitsche parameter gamma0

    public:

        /**
         * @brief Constructor that initializes the Nitsche parameters.
         * @param theta Nitsche parameter theta
         * @param gamma0 Nitsche parameter gamma0
         */
        NitscheContactEnforcement(scalar_type theta, scalar_type gamma0)
        : M_theta(theta), M_gamma0(gamma0) {}

        void enforce(getfem::model& md, const getfem::mesh_im& im) const override;
    };


    /**
     * @brief Penalty contact enforcement strategy.
     * This class implements the penalty method for enforcing contact conditions.
     * The penalty method applies a penalty term to the contact conditions.
     */
    class PenaltyContactEnforcement : public ContactEnforcementStrategy {

        scalar_type M_epsilon; ///< Penalty parameter epsilon

    public:

        /**
         * @brief Constructor for the PenaltyContactEnforcement class.
         * @param epsilon Penalty parameter epsilon
         */
        PenaltyContactEnforcement(scalar_type epsilon)
        : M_epsilon(epsilon) {}
    

        void enforce(getfem::model& md, const getfem::mesh_im& im) const override;
    };


} // namespace gf



#endif // _CONTACT_ENFORCEMENT_STRATEGY_HPP_