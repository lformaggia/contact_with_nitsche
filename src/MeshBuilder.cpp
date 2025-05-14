#include "MeshBuilder.hpp"
#include "MeshRegion.hpp"

bool DEBUGIR = false;
namespace gf {

    void BuiltInBuilder::buildMesh(getfem::mesh& mesh) const
    {
        size_type N = M_datafile("domain/problemDimension", 3);
        size_type Lx = M_datafile("domain/Lx", 2); // lenghts
        size_type Ly = M_datafile("domain/Ly", 4);
        size_type Lz = M_datafile("domain/Lz", 4);
        size_type Nx = M_datafile("domain/Nx", 4); // number of subdivisions
        size_type Ny = M_datafile("domain/Ny", 16);
        size_type Nz = M_datafile("domain/Nz", 16);
        std::string meshType = M_datafile("domain/meshType", "GT_PK(3,1)"); // element type for meshing, linear by defaults

        auto pgt = bgeot::geometric_trans_descriptor(meshType);

        if (true){
            auto pbDim = pgt->dim();
            assert(N==pbDim);
        }

        std::vector<size_type> nsubdiv(N);
        std::vector<scalar_type> lengths(N);
        nsubdiv[0] = Nx;
        nsubdiv[1] = Ny;
        nsubdiv[2] = Nz;
        lengths[0] = Lx;
        lengths[1] = Ly;
        lengths[2] = Lz;
    
        // create the mesh
        getfem::regular_unit_mesh(mesh, nsubdiv, pgt);

        // transform the mesh (scale the unit mesh to [Lx, Ly, Lz])
        //!\todo: maybe I also need to translate it, I want the fault to be positioned in x=0
        bgeot::base_matrix M(N,N);
        for (size_type i = 0; i < N; ++i)
            M(i,i) = lengths[i];

        mesh.transformation(M);
        std::cout << "Lx = " << Lx << ", -Lx = " << -Lx << std::endl;
        base_small_vector v {-0.5*Lx,0.,0.};
        mesh.translation(v);

        //!\todo: need to be sure that i can build the fracture according to the normal ...

    }


    void BuiltInBuilder::initRegions(const getfem::mesh& mesh, RegionMapType& regions) const {
        
        dal::bit_vector leftConvexesList;
        dal::bit_vector rightConvexesList;
        dal::bit_vector centeredConvexesList;

        
        // loop over all convexes in mesh
        for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
            
            bool insertedLeft = false, insertedRight = false, insertedCentered = false;

            // get the coords of the points of convex i
            auto local_points = mesh.points_of_convex(i);
            
            if (DEBUGIR){
                // std::clog << "Element: " << i << std::endl;
                // std::clog << "--Points\n";
            }

            // loop over the points of a convex
            for (const auto &pt: local_points){
                if (DEBUGIR)
                    std::clog <<"["<<pt[0]<<","<<pt[1]<<","<<pt[2]<<"]\n";
                
                if (pt[0] < -1e-6 && !insertedLeft){ // if one point has x < 0 --> left bulk region
                    leftConvexesList.add(i);
                    insertedLeft = true;
                }
                else if (pt[0] > 1.e-6 && !insertedRight){ // "   "   "   x > 0 --> right   "   "
                    rightConvexesList.add(i);
                    insertedRight = true;
                }
                if (std::abs(pt[0]) < 1.e-6 && !insertedCentered){
                    if (DEBUGIR) {
                        std::clog << "Adding an element in Fault convexList: "
                        << "pt[0] = " << pt[0]
                        << ", abs(pt[0]) = " << std::abs(pt[0]) << " < 1.e-6? "
                        << std::boolalpha << (std::abs(pt[0]) < 1.e-6) << std::endl;
                    }
                    centeredConvexesList.add(i);
                    insertedCentered = true;
                }
            }
        }
        regions["BulkLeft"] = std::make_unique<Bulk>(mesh, leftConvexesList, SideType::LEFT);
        regions["BulkRight"] = std::make_unique<Bulk>(mesh, rightConvexesList, SideType::RIGHT);
        regions["Fault"] = std::make_unique<Fault>(mesh,centeredConvexesList);

    }


    // GmshBuilder:: IMPLEMENTATION 

    void GmshBuilder::buildMesh(getfem::mesh& mesh) const
    {
        mesh.read_from_file(M_meshFile);
    }


    void GmshBuilder::initRegions(const getfem::mesh& mesh, RegionMapType &regions) const {
        std::map<std::string, size_type> regionMap;
        // getfem::import_mesh_gmsh(M_meshFile, mesh, regionMap);

        /* ... */
    
        auto itLeft = regionMap.find("BulkLeft");
        auto itRight = regionMap.find("BulkRight");
        auto itFault = regionMap.find("Fault");
    
        if (itLeft != regionMap.end()) {
            auto view = std::make_unique<Bulk>(mesh, itLeft->second, SideType::LEFT);
            regions["BulkLeft"] = std::move(view);
        }
    
        if (itRight != regionMap.end()) {
            auto view = std::make_unique<Bulk>(mesh, itRight->second, SideType::RIGHT);
            regions["BulkRight"] = std::move(view);
        }
    
        if (itFault != regionMap.end()) {
            auto view = std::make_unique<Fault>(mesh, itFault->second);
            regions["Fault"] = std::move(view);
        }
    }

}