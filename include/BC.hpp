#include "Core.hpp"
#include "Mesh.hpp"

namespace gf{


    class BC {
    
    protected:
        getfem::mesh_region M_region; ///< the mesh region where to apply the BC
        VectorFunctionType M_function; ///< The function
        BCType M_BCtype; ///< Dirichlet, Neumann or Mixed
        size_type M_ID; ///< the ID of the boundary face where the BC is applied

    public:

        BC(const getfem::mesh_region& rg, size_type regionID, VectorFunctionType func, BCType bctype)
        : M_region(rg), M_function(func), M_BCtype(bctype), M_ID(regionID) {}

        /**
         * @brief Returns the region (read only)
         */
        const getfem::mesh_region& getRegion() const { return M_region; };

        /**
         * @brief Return the BC type
         */
        virtual BCType type() const = 0;

        virtual std::string name() const = 0;

        /**
         * @brief Returns the ID of the face where the region is applied,
         * accordingly to the map built from BCHandler
         */
        size_type ID() const { return M_ID; }

        virtual ~BC() = default;

        /**
         * @brief Evaluate the function
         * @param x The node
         * @param t The time instant
         */
        base_small_vector eval(const base_node& x, scalar_type t) const {
            return M_function(x,t);
        }

        VectorFunctionType& f() { return M_function; }
        
        bool isLeftBoundary() const {
            return M_ID == 1 || M_ID == 2 || M_ID == 3 || M_ID == 4 || M_ID == 9;
        }
        
    };




    class BCDir : public BC {

    public:
        BCDir(const getfem::mesh_region& rg, size_type ID, VectorFunctionType f, BCType bctype)
        : BC(rg, ID, f, bctype){}

        BCType type() const override { return BCType::Dirichlet; }

        std::string name() const override {return "DirichletData"+std::to_string(M_ID); }

    };

    class BCNeu : public BC {

        bool M_isNormal;

    public:
        BCNeu(const getfem::mesh_region& rg, size_type ID, VectorFunctionType f, BCType bctype)
        : BC(rg, ID, f, bctype){}

        BCType type() const override { return BCType::Neumann; }

        std::string name() const override { return "NeumannData"+ std::to_string(M_ID);}

        bool isNormal() const { return M_isNormal; }

    };

} // namespace gf