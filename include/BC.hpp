#include "Core.hpp"
#include "Mesh.hpp"

namespace gf{


    class BC {
    
    protected:
        getfem::mesh_region M_region; ///< the mesh region where to apply the BC
        VectorFunctionType M_function; ///< The function
        BCType M_BCtype; ///< Dirichlet, Neumann or Mixed
        size_type M_ID; ///< the ID of the boundary face where the BC is applied
        base_small_vector M_normal;

    public:

        BC(const getfem::mesh_region& rg, size_type regionID,
            VectorFunctionType func, BCType bctype, base_small_vector n)
        : M_region(rg), M_function(func), M_BCtype(bctype), M_ID(regionID), M_normal(n){
        }

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

        base_small_vector normal() const { return M_normal; }
        
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
        BCDir(const getfem::mesh_region& rg, size_type ID, VectorFunctionType f, BCType bctype, base_small_vector n)
        : BC(rg, ID, f, bctype, n){
        }

        BCType type() const override { return BCType::Dirichlet; }

        std::string name() const override {return "DirichletData"+std::to_string(M_ID); }

    };

    class BCNeu : public BC {

    public:
        BCNeu(const getfem::mesh_region& rg, size_type ID, VectorFunctionType f, BCType bctype, base_small_vector n)
        : BC(rg, ID, f, bctype, n){}

        BCType type() const override { return BCType::Neumann; }

        std::string name() const override { return "NeumannData"+ std::to_string(M_ID);}

    };


    class BCMix : public BC {
    public:
        BCMix(const getfem::mesh_region& rg, size_type ID, VectorFunctionType f, BCType bctype, base_small_vector n)
        : BC(rg, ID, f, bctype, n){}

        BCType type() const override { return BCType::Mixed; }

        std::string name() const override { return "MixedData"+ std::to_string(M_ID);}

    };

} // namespace gf