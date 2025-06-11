#include "Mesh.hpp"

bool DEBUG = false;

namespace gf {

    Mesh::Mesh(const Params& p)
    : M_params(p)
    {
        if (DEBUG)
            std::clog << "============= MESH INFORMATION =============\n";

        if (M_params.gmsh)
            M_meshBuilder = std::make_unique<GmshBuilder>(M_params.domain);
        else
            M_meshBuilder = std::make_unique<BuiltInBuilder>(M_params.domain);
        
        M_meshBuilder->construct(M_mesh);
            
        if (DEBUG){
            std::clog << "\n\n============= REGION INFORMATION =============" << std::endl;
            M_mesh.write_to_file(std::cout);

            // Check boundary regions
            {
                std::clog << std::endl;
                for (size_type i = 1; i < 11; ++i) {
                    const getfem::mesh_region& region = M_mesh.region(i);

                    std::clog << "Boundary region " << i << " has " << region.size() << " elements" << std::endl;
                    std::clog << "Number of convexes in region: " << region.nb_convex() << std::endl;
                    std::clog << "Region contains only faces: " << region.is_only_faces() << std::endl;
                    std::clog << "Region contains only convexes: " << region.is_only_convexes() << std::endl;

                    std::clog << "Elements in region " << i << ": ";
                    for (getfem::mr_visitor it(region); !it.finished(); ++it) {
                        std::clog << "(" << it.cv() << ", " << it.f() << ") ";
                    }
                    std::clog << std::endl << std::endl;
                }
            }

            // Check BulkLeft
            {
                const auto& region = M_mesh.region(BulkLeft);
                std::clog << "BULKLEFT REGION has " << region.size() << " elements" << std::endl;
                std::clog << "Number of convexes in region: " << region.nb_convex() << std::endl;
                std::clog << "Region contains only faces: " << region.is_only_faces() << std::endl;
                std::clog << "Region contains only convexes: " << region.is_only_convexes() << std::endl;

                std::clog << "Elements in region BulkLeft: ";
                for (getfem::mr_visitor it(region); !it.finished(); ++it) {
                    std::clog << it.cv() << " ";
                }
                std::clog << std::endl << std::endl;
            }

            // Check BulkRight
            {
                const auto& region = M_mesh.region(BulkRight);
                std::clog << "BULKRIGHT REGION has " << region.size() << " elements" << std::endl;
                std::clog << "Number of convexes in region: " << region.nb_convex() << std::endl;
                std::clog << "Region contains only faces: " << region.is_only_faces() << std::endl;
                std::clog << "Region contains only convexes: " << region.is_only_convexes() << std::endl;

                std::clog << "Elements in region BulkRight: ";
                for (getfem::mr_visitor it(region); !it.finished(); ++it) {
                    std::clog << it.cv() << " ";
                }
                std::clog << std::endl << std::endl;
            }

            // Check Fault Region
            {
                const getfem::mesh_region &fault_region = M_mesh.region(Fault);
                std::clog << "FAULT REGION has " << fault_region.size() << " elements" << std::endl;
                std::clog << "Number of convexes in region: " << fault_region.nb_convex() << std::endl;
                std::clog << "Region contains only faces: " << fault_region.is_only_faces() << std::endl;
                std::clog << "Region contains only convexes: " << fault_region.is_only_convexes() << std::endl;

                // Loop through (convex, face) pairs
                std::clog << "Faces in region Fault (convex, face): ";
                for (getfem::mr_visitor it(fault_region); !it.finished(); ++it) {
                    std::clog << "(" << it.cv() << ", " << it.f() << ") ";
                }
                std::clog << std::endl;
            }            
        }

        
    }


} // namespace gf