#ifndef _FEM_MANAGER_HPP_
#define _FEM_MANAGER_HPP_

#include "Core.hpp"
#include "Params.hpp"

namespace gf {

    class MeshBuilderStrategy;
    
    /** \todo: Check if it's better to:
     * - do a forward declaration like this
     * - forward declare in Core.hpp --> define RegionMapType there
     * - include MeshBuilderStrategy here
     */

    class FEMManager {
    public:
        FEMManager(const getfem::mesh&);
    
        void setMeshFem(const Numerics& n, const getfem::mesh&);

        getfem::mesh_fem& mf_u() { return M_mfU; }
        
        getfem::mesh_fem& mf_stress() { return M_mfStress; }

        getfem::mesh_fem& mf_rhs() { return M_mfRhs; }

        getfem::mesh_fem& mf_LM() { return M_mfLM; }

    private:
        getfem::mesh_fem M_mfU;
        getfem::mesh_fem M_mfStress;
        getfem::mesh_fem M_mfRhs;
        getfem::mesh_fem M_mfLM;

    };

}



#endif // _FEM_MANAGER_HPP_