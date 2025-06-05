#include "GetPot"
#include "Mesh.hpp"
#include "ContactProblem.hpp"

int main(int argc, char * argv[]){

    GETFEM_MPI_INIT(argc, argv);

    using namespace gf;

    // parse command line options
    GetPot command_line(argc, argv);
    const std::string dataFileName = command_line.follow("data.pot", 2, "-f",
            "--file");
    const std::string meshFileName = command_line.follow("", 2, "-m",
            "--mesh");
    const bool verbose = command_line.search("-v");
    
    std::ifstream datafile(dataFileName);
    if(!datafile.is_open())
        throw std::runtime_error("Could not open the datafile!");
    
    std::ifstream meshfile(meshFileName);
    if(!meshfile.is_open())
        std::cerr << "Could not open the mesh file! Using the datafile parameters..." << std::endl;

    Params p(dataFileName, meshFileName, verbose);

    Mesh mesh(p);
        
    ContactProblem pb(mesh, p);
    
    pb.init();

    pb.assemble();

    pb.solve();

    GETFEM_MPI_FINALIZE;

    return 0;

}