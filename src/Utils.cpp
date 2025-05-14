#include "Utils.hpp"

namespace gf {

    // Remove brackets
    std::vector<std::string>
    splitString(std::string stringValue) {
        stringValue.erase(std::remove(stringValue.begin(), stringValue.end(), '['), stringValue.end());
        stringValue.erase(std::remove(stringValue.begin(), stringValue.end(), ']'), stringValue.end());

        // Split by commas
        std::vector<std::string> components;
        std::istringstream ss(stringValue);
        std::string item;

        while (std::getline(ss, item, ',')) {
            components.push_back(item);
        }

        return components;
    }

    std::vector<size_type>
    toVec(std::string stringValue) {
        auto components = splitString(stringValue);
        std::vector<size_t> result; 
        for (const auto& c:components)
            result.emplace_back(std::stod(c));
        return result;
    }

} // namespace gf
