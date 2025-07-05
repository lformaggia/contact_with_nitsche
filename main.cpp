#include "Mesh.hpp"
#include "ContactProblem.hpp"
#include "Utils.hpp"

int main(int argc, char * argv[]){
    using namespace gf;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        }
    }

    Params p(argc, argv);
    Mesh mesh(p);
    ContactProblem pb(mesh, p);
    
    pb.init();
    pb.assemble();
    pb.solve();

    return 0;

}