#ifndef _CORE_HPP_
#define _CORE_HPP_

// Level Set and Xfem stuff:
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_config.h>
#include <getfem/getfem_assembling.h> // import assembly methods (and comp. of 
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>     // export functions (save the solution in a file)
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_MUMPS_interface.h>
#include <gmm/gmm_superlu_interface.h>

#include <getfem/bgeot_mesh.h>

#include <string>

#include <getfem/getfem_mesh_fem_product.h>
#include <getfem/getfem_mesh_fem_global_function.h>

#include "GetPot"
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <iomanip>
#include <functional>


/* some Getfem++ types that we will be using */
namespace gf{
    using bgeot::scalar_type; ///< = double
    using bgeot::base_small_vector; ///< special class for small (dim < 16) vectors
    using bgeot::base_node; ///< geometrical nodes (derived from base_small_vector)
    using bgeot::size_type; ///< basically std::vector<scalar_type>::size_type

    using ScalarFunctionType = std::function<scalar_type(base_node,scalar_type)> ///< f(\vector{x},t)
    using VectorFunctionType = std::function<base_small_vector(base_node, scalar_type)> ///< \vector{f}(\vector{x},t)

    enum SideType{
        LEFT,
        RIGHT
    };

    enum BCType {
        Dirichlet,
        Neumann,
        Mixed
    }


    /* definition of some matrix/vector types. These ones are built
    using the predefined types in Gmm++ */

    // using sparseVector_Type = gmm::rsvector<scalar_type>;  // sparse std::vector optimized for read operations

    // typedef gmm::row_matrix<sparseVector_Type> sparseMatrix_Type;  // row_matrix
    // typedef boost::shared_ptr<sparseMatrix_Type> sparseMatrixPtr_Type;  // ptr to a row-matrix
    // typedef std::vector<sparseMatrix_Type> sparseMatrixContainer_Type;  // vector of row_matrices
    // typedef std::vector<sparseMatrixPtr_Type> sparseMatrixPtrContainer_Type;  // vector of ptrs to row-matrix

    // typedef std::vector<scalar_type> scalarVector_Type;  // std::vector
    // typedef boost::shared_ptr<scalarVector_Type> scalarVectorPtr_Type;  // ptr to std::vector
    // typedef std::vector<scalarVector_Type> scalarVectorContainer_Type;  // std::vector of std::vectors
    // typedef boost::shared_ptr<scalarVectorPtr_Type> scalarVectorContainerPtr_Type; // ptr to ptr to std::vector
    // typedef std::vector<scalarVectorPtr_Type> scalarVectorPtrContainer_Type;  // vectir of ptrs to vector

    // typedef std::vector<size_type> sizeVector_Type;  // size_type
    // typedef boost::shared_ptr<sizeVector_Type> sizeVectorPtr_Type;  // ptr to size_type
    // typedef std::vector<sizeVector_Type> sizeVectorContainer_Type;  // vector of size_type
    // typedef std::vector<sizeVectorPtr_Type> sizeVectorPtrContainer_Type;  // vector of ptr

    // typedef std::pair<size_type, size_type> pairSize_Type;
    // typedef std::vector<pairSize_Type > pairSizeVector_Type;
    // typedef std::vector<pairSizeVector_Type> pairSizeVectorContainer_Type;

    // typedef std::vector < std::string > stringContainer_Type;

    // typedef boost::shared_ptr<getfem::mesh> GFMeshFEMPtr_Type;
    // typedef std::vector<GFMeshFEMPtr_Type> GFMeshFEMPtrContainer_Type;

    // namespace LifeV{
    // typedef scalar_type Double;
    // typedef size_type UInt;
    // typedef size_type ID;
    // typedef int Int;
    // }
    // enum ElementDimension
    // {
    //     MEDIUM = 3, FRACTURE = MEDIUM-2
    // };

}
#endif // _CORE_HPP_
