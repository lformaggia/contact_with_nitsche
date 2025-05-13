#include "Core.hpp"

namespace gf{


    class BC {
    
    protected:
        getfem::mesh_region M_region; ///< the mesh region where to apply the BC
        VectorFunctionType M_function; ///< The function
        BCType M_BCtype;
        size_type M_ID;

    public:

        BC(const getfem::mesh_region&, size_type, VectorFunctionType&&, BCType);

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
        BCDir(const getfem::mesh_region& region, size_type ID, VectorFunctionType&& f, BCType bctype)
        : BC(region, ID, std::move(f), bctype){}

        std::string type() const override { return "Dirichlet"; }

    };

    class BCNeu : public BC {

        bool M_isNormal;

    public:
        BCNeu(const getfem::mesh_region& region, size_type ID, VectorFunctionType&& f, BCType bctype)
        : BC(region, ID, std::move(f), bctype){}

        std::string type() const override { return "Neumann"; }

        bool isNormal() const { return M_isNormal; }

    };

} // namespace gf