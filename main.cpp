#include "Mesh.hpp"
#include "ContactProblem.hpp"

int main(int argc, char * argv[]){

    using namespace gf;

    Params p(argc, argv);
    Mesh mesh(p);
    ContactProblem pb(mesh, p);
    
    pb.init();
    pb.assemble();
    pb.solve();

    return 0;

}