#include "GetPot"
#include "ContactProblem.hpp"

int main(int argc, char * argv[]){
    using namespace gf;

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
        
    ContactProblem pb(dataFileName, meshFileName, verbose);
    
    pb.init();



    

    return 0;

}