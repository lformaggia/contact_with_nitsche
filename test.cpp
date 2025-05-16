#include "GetPot"
#include <string>
#include <sstream>
#include <vector>

std::vector<std::size_t>
toVec(std::string stringValue);


#include "mpParser.h"
#include "muParserXInterface.hpp"
#include <iostream>
#include <functional>
#include <vector>
#include <memory>




int main(int argc, char * argv[]){

    // parse command line options
    GetPot command_line(argc, argv);
    const std::string dataFileName = command_line.follow("data.pot", 2, "-f",
            "--file");
    const std::string meshFileName = command_line.follow("", 2, "-m",
            "--mesh");
    const bool verbose = command_line.search("-v");
    
    std::ifstream test(dataFileName);
    if(!test.is_open())
        throw std::runtime_error("Could not open the file!");
    
    
    GetPot datafile(dataFileName.c_str());

    std::string regionDisp = datafile("physics/regionDisp", "");
    std::vector<std::size_t> regionsID = toVec(regionDisp);
    
    std::cout << "TEST\nRegionsID: ";
    for (const auto&v: regionsID) std::cout << v << " ";
    std::cout << std::endl;

    using namespace MuParserInterface2;
    muParserXInterface mx;

    std::vector<double> x = {1, 2, 3};
    double t = 100;
    std::vector<std::string> setexpr{"x[0]+sin(x[1])", "0.", "x[1]^2"};
    // for (auto& expr: setexpr){
    //     mx.set_expression(expr);
    //     std::cout << "Evaluating: mx(x) = ("
    //     << mx(x)[0] << std::endl;
    //     expr.clear();
    // }
    std::string expr = "{ x[0] + sin(x[1]), sqrt(t)+x[2], x[1]^2 }";
    mx.set_expression(expr);
    std::cout << "Evaluating: mx(x) = ("
        << mx(x,t)[0] << ", " << mx(x,t)[1] << ", " << mx(x,t)[2] << ")" << std::endl;
    expr.clear();

    using VectorFunctionType = std::function<std::vector<double>(std::vector<double>, double)>;
    
    expr = "{ x[0] + log(x[2]), 0., 3*x[1] + t}";
    mx.set_expression(expr);
    VectorFunctionType myFunc(mx);
    std::cout << "Evaluating: myFunc(x) = ("
        << myFunc(x,t)[0] << ", " << myFunc(x,t)[1] << ", " << myFunc(x,t)[2] << ")" << std::endl;
    expr.clear();

    // expr = "{ x[0] , x[1], x[2]}";
    // mx.set_expression(expr);
    // VectorFunctionType myFunc2(mx);
    // std::cout << "Evaluating: myFunc2(x) = ("
    //     << myFunc2(x)[0] << ", " << myFunc2(x)[1] << ", " << myFunc2(x)[2] << ")" << std::endl;
    // expr.clear();

    return 0;

}

std::vector<std::size_t>
toVec(std::string stringValue) {
    stringValue.erase(std::remove(stringValue.begin(), stringValue.end(), '['), stringValue.end());
    stringValue.erase(std::remove(stringValue.begin(), stringValue.end(), ']'), stringValue.end());

    // Split by commas
    std::vector<std::string> components;
    std::istringstream ss(stringValue);
    std::string item;

    while (std::getline(ss, item, ',')) {
        components.push_back(item);
    }
    std::vector<size_t> result; 
    for (const auto& c:components)
        result.emplace_back(std::stod(c));

    return result;
}