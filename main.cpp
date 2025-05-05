#include "GetPot"
#include "ContactProblem.hpp"

int main(int argc, char * argv[]){
    using namespace gf;

    GetPot command_line(argc, argv);
    const std::string data_file_name = command_line.follow("data", 2, "-f",
            "--file");

    GetPot dataFile(data_file_name.data());
    Domain myDom(dataFile);

    return 0;

}