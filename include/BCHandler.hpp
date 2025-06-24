#ifndef _BC_HANDLER_HPP_
#define _BC_HANDLER_HPP_

#include "Core.hpp"
#include "Params.hpp"
#include "Mesh.hpp"
#include "GetPot"
#include "BC.hpp"
#include "muParserXInterface.hpp"

#include <unordered_map>

namespace gf {

    /**
     * @brief The BCHandler class is responsible for reading and managing boundary conditions
     * from a configuration file, and providing access to the different types of boundary conditions.
     */
    class BCHandler {

        using BCListType = std::unordered_map<BCType, std::vector<std::unique_ptr<BC>>>;
        using BCstringsType = std::unordered_map<BCType, std::vector<std::string>>;

        const getfem::mesh& M_mesh; ///< The mesh object
        BCListType M_BCList; ///< The BCobjects
        BCstringsType M_BCStrings; ///< The strings of the prescribed BC functions (just for output)
        
        muParserXInterface M_parser; ///< The muParserX interface for parsing expressions

        /**
         * @brief Reads a boundary condition of type T from the GetPot configuration.
         * @tparam T The type of boundary condition to read (Dirichlet, Neumann, Mixed).
         * @param pot The GetPot object containing the BC data.
         */
        template <BCType T>
        void read(const GetPot&);

    public:

        /**
         * @brief Constructor that initializes the BCHandler with a mesh.
         * @param m The mesh object to associate with this BCHandler.
         */
        BCHandler(const getfem::mesh& m);

        /**
         * @brief Reads boundary conditions from a GetPot configuration file.
         * @param pot The GetPot object containing the configuration data.
         */
        void readBC(const GetPot& pot);

        /**
         * @brief Getter for the Neumann boundary conditions.
         * @return A constant reference to a vector of unique pointers to Neumann boundary conditions.
         * @note The vector contains unique pointers to BC objects, which are polymorphic.
         */
        const std::vector<std::unique_ptr<BC>> & Neumann() const;

        /**
         * @brief Getter for the Dirichlet boundary conditions.
         * @return A constant reference to a vector of unique pointers to Dirichlet boundary conditions.
         * @note The vector contains unique pointers to BC objects, which are polymorphic.
         */
        const std::vector<std::unique_ptr<BC>> & Dirichlet() const;

        /**
         * @brief Getter for the Mixed boundary conditions.
         * @return A constant reference to a vector of unique pointers to Mixed boundary conditions.
         * @note The vector contains unique pointers to BC objects, which are polymorphic.
         */
        const std::vector<std::unique_ptr<BC>> & Mixed() const;

    };

} // namespace gf



#endif // _BC_HANDLER_HPP_