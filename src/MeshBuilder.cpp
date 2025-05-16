#include "MeshBuilder.hpp"
#include "MeshRegion.hpp"

bool DEBUGIR = false;
namespace gf {

    void BuiltInBuilder::buildMesh(getfem::mesh& mesh) const
    {
        size_type N = M_domain.dim;
        size_type Nx = M_domain.Nx; // number of subdivisions
        size_type Ny = M_domain.Ny;
        size_type Nz = M_domain.Nz;
        std::string meshType = M_domain.meshType; // element type for meshing, linear by defaults

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
        lengths[0] = M_domain.Lx;
        lengths[1] = M_domain.Ly;
        lengths[2] = M_domain.Lz;
    
        // create the mesh
        getfem::regular_unit_mesh(mesh, nsubdiv, pgt);

        // transform the mesh (scale the unit mesh to [M_domain.Lx, M_domain.Ly, M_domain.Lz])
        //!\todo: maybe I also need to translate it, I want the fault to be positioned in x=0
        bgeot::base_matrix M(N,N);
        for (size_type i = 0; i < N; ++i)
            M(i,i) = lengths[i];
        mesh.transformation(M);

        base_small_vector v {-0.5*M_domain.Lx,0.,0.};
        mesh.translation(v);

        //!\todo: need to be sure that i can build the fracture according to the normal ...

    }


    void
    BuiltInBuilder::initRegions(const getfem::mesh& mesh,
                                RegionMapType& regions,
                                BoundaryMapType& bds) const
        {
        // Create the internal regions BulkLeft, BulkRight, Fault
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

        
        // Create the boundary regions (map Id -> boundary region)
        getfem::mesh_region borderFaces;
        getfem::outer_faces_of_mesh(mesh, borderFaces);

        bds[1] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{-M_domain.Lx/2.,0.,0.},
            base_node{0.,M_domain.Ly,0.}),
            1
        );
        bds[2] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{0.,M_domain.Ly,0.},
            base_node{-M_domain.Lx/2,M_domain.Ly,M_domain.Lz}),
            2
        );
        bds[3] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{-M_domain.Lx/2,M_domain.Ly,M_domain.Lz},
            base_node{0.,0.,M_domain.Lz}),
            3
        );
        bds[4] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{0.,0.,M_domain.Lz},
            base_node{-M_domain.Lx/2.,0.,0.}),
            4
        );
        bds[5] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{0.,0.,0.},
            base_node{M_domain.Lx/2.,M_domain.Ly,0.}),
            5
        );
        bds[6] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{M_domain.Lx/2.,M_domain.Ly,0.},
            base_node{0.,M_domain.Ly,M_domain.Lz}),
            6
        );
        bds[7] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{0.,M_domain.Ly,M_domain.Lz},
            base_node{M_domain.Lx/2.,0.,M_domain.Lz}),
            7
        );
        bds[8] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{M_domain.Lx/2.,0.,M_domain.Lz},
            base_node{0.,0.,0.}),
            8
        );
        
        // alternatively:
        // bds(9) = getfem::select_faces_of_normal(mesh, border_faces,
        //    base_small_vector(-1.,0.,0.), 0.01);
        // bds(10) = getfem::select_faces_of_normal(mesh, border_faces,
        //    base_small_vector(1.,0.,0.), 0.01);
        bds[9] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{-M_domain.Lx/2.,0.,0.},
            base_node{-M_domain.Lx/2.,M_domain.Ly,M_domain.Lz}),
            9
        );
        bds[10] = std::make_unique<Boundary> (mesh, getfem::select_faces_in_box(
            mesh, borderFaces,
            base_node{M_domain.Lx/2.,0.,0},
            base_node{M_domain.Lx/2.,M_domain.Ly,M_domain.Lz}),
            10
        );

        // ALTERNATIVE:
        // for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
        //     assert(it.is_face());
        //     base_node x = /* get the centroid of the face*/

        //     base_small_vector un = mesh.normal_of_face_of_convex(it.cv(), it.f());
        //     un /= gmm::vect_norm2(un);

        //     if ((un[2] + 1.0) < 1.0E-7)
        //         M_IdToRegion[1].add(it.cv(), it.f());
        //     else if ((un[1] - 1.0 < 1.0E-7) && x[0] < 0) 
        //         M_IdToRegion[2].add(it.cv(),it.f());
        //     else if ((un[2] - 1.0 < 1.0E-7) && x[0] < 0)
        //         M_IdToRegion[3].add(it.cv(),it.f());
        //     else if ((un[1] + 1.0 < 1.0E-7) && x[0] < 0)
        //         M_IdToRegion[4].add(it.cv(),it.f());
        //     else if ((un[2] + 1.0) < 1.0E-7)
        //         M_IdToRegion[5].add(it.cv(), it.f());
        //     else if ((un[1] - 1.0 < 1.0E-7) && x[0] > 0) 
        //         M_IdToRegion[6].add(it.cv(),it.f());
        //     else if ((un[2] - 1.0 < 1.0E-7) && x[0] > 0)
        //         M_IdToRegion[7].add(it.cv(),it.f());
        //     else if ((un[1] + 1.0 < 1.0E-7) && x[0] > 0)
        //         M_IdToRegion[8].add(it.cv(),it.f());
        //     else if (un[0] + 1.0 < 1.0E-7)
        //         M_IdToRegion[9].add(it.cv(),it.f());
        //     else if (un[0] - 1.0 < 1.0E-7)
        //         M_IdToRegion[10].add(it.cv(),it.f());
        // }

    }


    // GmshBuilder:: IMPLEMENTATION 
    /** \todo */
    void
    GmshBuilder::buildMesh(getfem::mesh& mesh) const
    {
        mesh.read_from_file(M_meshFile);
    }


    /** \todo */
    void
    GmshBuilder::initRegions(const getfem::mesh& mesh,
                            RegionMapType &regions,
                            BoundaryMapType &bds) const
    {
        std::map<std::string, size_type> regionMap;
        // getfem::import_mesh_gmsh(M_meshFile, mesh, regionMap);

        /** \todo ... */
    
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