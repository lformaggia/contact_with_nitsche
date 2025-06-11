#include "Mesh.hpp"
#include "ContactProblem.hpp"

int main(int argc, char * argv[]){

    // GETFEM_MPI_INIT(argc, argv);

    using namespace gf;

    Params p(argc, argv);
    Mesh mesh(p);
    ContactProblem pb(mesh, p);
    
    pb.init();
    pb.assemble();
    pb.solve();

    // GETFEM_MPI_FINALIZE;

    return 0;

}