#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include "Core.hpp"

namespace gf {

    /**
     * @brief Converts a string of type [1,7,9,...] to a vector of strings.
     * The string should contain comma-separated values.
     * Used inside toVec().
     * @param stringValue The input string to be converted.
     * @return A vector containing the converted values.
     */
    std::vector<std::string> splitString(std::string stringValue);

    /**
     * @brief Converts a string of type [1,7,9,...] to a vector of size_type.
     * The string should contain comma-separated values.
     * Used to read the region IDs for applying boundary conditions from the datafile.
     * @param stringValue The input string to be converted.
     * @return A vector containing the converted values as size_type.
     */
    std::vector<size_type> toVec(std::string stringValue);

    /**
     * @brief Print optional runtime options when running
     *      ./main -h
     */
    void print_help();

} // namespace gf

#endif // _UTILS_HPP_