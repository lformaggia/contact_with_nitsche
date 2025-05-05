#include "Core.hpp"

#include <functional>

namespace gf{


    class BC {
    
    protected:
        getfem::mesh & M_mesh;
        getfem::mesh_region M_border_faces;
        getfem::mesh_region M_region; ///< the mesh region where to apply the BC
        VectorFunctionType M_function; ///< The function
        std::string M_expr; ///< The expression (GWFL)
        BCType M_BCtype;

    public:
        /** 
         * @brief Constructs the BC object (Neumann, Dirichlet)
         * @param datafile The datafile
         */ 
        BC(const GetPot& datafile);

        /**
         * @brief Returns the region (read only)
         */
        const getfem::mesh_region& getRegion() const { return M_region; };

        /**
         * @brief Return the BC type
         */
        virtual BCType type() const = 0;

        virtual ~BC() = default;

        /**
         * @brief Evaluate the function
         * @param x The node
         * @param t The time instant
         */
        small_base_vector eval(const base_node& x, scalar_type t) const;
        
    }



    class BCDir : public BC {

    public:
        std::string type() const override { return "Dirichlet"; }

    };

    class BCNeu : public BC {
    private:
        bool isNormal;
    public:
        BCneu(const GetPot& datafile);

        std::string type() const override { return "Neumann"; }

        std::string isNormal() const { return isNormal; }

    };

} // namespace gf