#ifndef _CORE_HPP_
#define _CORE_HPP_


#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_config.h>
#include <getfem/getfem_assembling.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_MUMPS_interface.h>
#include <gmm/gmm_superlu_interface.h>

#include <getfem/bgeot_mesh.h>

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <memory>

#include <getfem/getfem_mesh_fem_product.h>
#include <getfem/getfem_mesh_fem_global_function.h>
#include <getfem/getfem_models.h>
#include <getfem/getfem_model_solvers.h>
#include <getfem/getfem_plasticity.h>

#include "GetPot"
#include <iostream>
#include <iomanip>
#include <functional>

namespace gf {
    
    // Some Getfem++ types that are used
    using bgeot::scalar_type; ///< = double
    using bgeot::base_small_vector; ///< special class for small (dim < 16) vectors
    using bgeot::base_vector; ///< just a std::vector<scalar_type>
    using bgeot::base_node; ///< geometrical nodes (derived from base_small_vector)
    using bgeot::size_type; ///< basically std::vector<scalar_type>::size_type
    using bgeot::dim_type; ///< type for the dimension of the mesh

    using ScalarFunctionType = std::function<scalar_type(base_node,scalar_type)>; ///< f(\vector{x},t)
    using VectorFunctionType = std::function<base_small_vector(base_node, scalar_type)>; ///< \vector{f}(\vector{x},t)
    
    using plain_vector = getfem::modeling_standard_plain_vector; /// a vector of scalar_type, used for Getfem++ models
    using row_matrix = gmm::row_matrix<plain_vector>; 

    /**
     * @brief Enumeration for boundary condition types.
     * This enum defines the types of boundary conditions that can be applied to the model.
     */
    enum BCType {
        Dirichlet,
        Neumann,
        Mixed
    };

    /**
     * @brief Enumeration for regions in the mesh.
     * This enum defines the ID associated with different regions in the mesh.
     * These regions are used to apply boundary conditions or to identify specific parts of the mesh.
     */
    enum RegionType {
        /* This is the order in which the physical entities are written out on the .geo file.
            BottomLeft = 1,
            YmaxLeft = 2,
            TopLeft = 3,
            YminLeft = 4,
            BottomRight = 5,
            YmaxRight = 6,
            TopRight = 7,
            YminRight = 8,
            Xmin = 9,
            Xmax = 10,
        */
        BulkLeft = 11,
        BulkRight = 12,
        Fault = 13
    };

} // namespace gf

#endif // _CORE_HPP_

