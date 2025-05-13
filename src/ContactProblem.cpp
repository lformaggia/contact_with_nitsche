#include "ContactProblem.hpp"
bool DEBUG = true;
namespace gf{

    ContactProblem::ContactProblem (const std::string& filename, const std::string& meshfile, bool ver)
    : M_datafile(filename.c_str()), M_params(M_datafile, meshfile,ver){
    }

    void ContactProblem::init() {
        if (M_params.meshFile.empty()){
            M_meshBuilder = std::make_unique<BuiltInBuilder>(M_datafile);
            if (M_params.verbose)
                std::cout << "Building the mesh internally... ";
        }
        else {
            M_meshBuilder = std::make_unique<GmshBuilder>(M_params.meshFile);
            std::cout << "Reading the mesh from " << M_params.meshFile << "... ";
        }

        M_meshBuilder->buildMesh(M_mesh);
        if (M_params.verbose)
            std::cout << "done." << std::endl;

        ////////////////////////
        if (DEBUG){
            M_mesh.write_to_file(std::clog);
        }
        ////////////////////////

            
        if (M_params.verbose)
            std::cout << "\nCreating bulk and fault regions... ";
        M_meshBuilder->initRegions(M_mesh, M_regions);
        if (M_params.verbose)
            std::cout << "done." << std::endl;

        ////////////////////////
        if (DEBUG){
            /** \todo: output information related to the mesh regions */
            std::clog << "Region Information:" << std::endl;
            for (const auto& [name,r] : M_regions) {
                std::clog << "Region " << r->region().index() << " size: " << r->region().size() << std::endl;
                std::clog << "Name " << name << std::endl;
                std::clog << "Number of convexes in region: " << r->region().nb_convex() << std::endl;
                std::cout << "Region " << r->region().index() << " contains only faces: " << r->region().is_only_faces() << std::endl;
                std::cout << "Region " << r->region().index() << " contains only convexes: " << r->region().is_only_convexes() << std::endl;
                const dal::bit_vector& convexes = r->region().index();
                std::cout << "Convexes in region " << name << ": ";
                for (dal::bv_visitor i(convexes); !i.finished(); ++i) {
                    std::cout << i << " ";
                }
                std::cout << std::endl;
                // base_node Pmin, Pmax;
                // r->region().bounding_box(Pmin, Pmax);
                // std::clog << "Bounding box of region: [" << Pmin << ", " << Pmax << "]" << std::endl;
                // std::cout << std::endl;

            }

        }
        ////////////////////////


        M_BC = std::make_unique<BCHandler>(M_mesh, M_params.domain);
        M_BC->readBC(M_datafile);
        // M_FEM.setMeshFem(M_datafile, M_mesh, M_regions);
    }

    

    

}