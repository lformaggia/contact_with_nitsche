#include "Core.hpp"

namespace gf{


    class BC {
    
    protected:
        getfem::mesh_region M_region; ///< the mesh region where to apply the BC
        VectorFunctionType M_function; ///< The function
        BCType M_BCtype;

    public:

        BC(const getfem::mesh_region&, size_type, VectorFunctionType&&);

        /**
         * @brief Returns the region (read only)
         */
        const getfem::mesh_region& getRegion() const { return M_region; };

        /**
         * @brief Return the BC type
         */
        virtual std::string type() const = 0;

        virtual ~BC() = default;

        /**
         * @brief Evaluate the function
         * @param x The node
         * @param t The time instant
         */
        base_small_vector eval(const base_node& x, scalar_type t) const {
            return M_function(x,t);
        }
        
    };




    class BCDir : public BC {

    public:

        std::string type() const override { return "Dirichlet"; }

    };

    class BCNeu : public BC {

        bool isNormal;

    public:

        std::string type() const override { return "Neumann"; }

        bool isNormal() const { return isNormal; }

    };

} // namespace gf