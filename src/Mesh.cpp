#include "Mesh.hpp"

bool DEBUG = true;

namespace gf {

    Mesh::Mesh(const Params& p)
    : M_params(p)
    {
        if (DEBUG)
            std::clog << "============= MESH INFORMATION =============\n";

        if (M_params.meshFile.empty()){
            M_meshBuilder = std::make_unique<BuiltInBuilder>(M_params.domain);
            if (p.verbose)
                std::clog << "Building the mesh internally... ";
        }
        else {
            M_meshBuilder = std::make_unique<GmshBuilder>(M_params.meshFile);
            if (p.verbose)
                std::clog << "Reading the mesh from " << M_params.meshFile << "... ";
        }
        
        M_meshBuilder->buildMesh(M_mesh);

        if (M_params.verbose)
            std::clog << "done." << std::endl;
            
        ////////////////////////
        if (DEBUG){
            M_mesh.write_to_file(std::clog);
            std::clog << "\n\n============= REGION INFORMATION =============" << std::endl;
        }
        ////////////////////////

            
        if (M_params.verbose)
            std::clog << "Creating bulk, fault and boundary regions... ";

        M_meshBuilder->initRegions(M_mesh, M_regions, M_bdRegions);

        if (M_params.verbose)
            std::clog << "done." << std::endl;

            
        ////////////////////////
        if (DEBUG){
            // std::clog << "Fault is boundary of Bulkleft?"
            //     << M_regions["Fault"]->region().region_is_faces_of() << std::endl;
            /** Output information related to the mesh regions */
            std::clog << std::endl;
            for (const auto& [name,r] : M_regions) {
                std::clog << "Region " << r->region().index() << " size: " << r->region().size() << std::endl;
                std::clog << "Name " << name << std::endl;
                std::clog << "Number of convexes in region: " << r->region().nb_convex() << std::endl;
                std::clog << "Region " << r->region().index() << " contains only faces: " << r->region().is_only_faces() << std::endl;
                std::clog << "Region " << r->region().index() << " contains only convexes: " << r->region().is_only_convexes() << std::endl;
                const dal::bit_vector& convexes = r->region().index();
                std::clog << "Convexes in region " << name << ": ";
                for (dal::bv_visitor i(convexes); !i.finished(); ++i) {
                    std::clog << i << " ";
                }
                std::clog << std::endl;
                // base_node Pmin, Pmax;
                // r->region().bounding_box(Pmin, Pmax);
                // std::clog << "Bounding box of region: [" << Pmin << ", " << Pmax << "]" << std::endl;
                // std::cout << std::endl;
            }
        }

        
    }


} // namespace gf