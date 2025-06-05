#include "MeshBuilder.hpp"

bool DEBUGIR = true;
namespace gf {

    void
    MeshBuilderStrategy::construct(getfem::mesh& mesh) const
    {
        std::cout << "Building the mesh internally... ";
        buildMesh(mesh);
        std::cout << "done.\n";
        std::cout << "Initializing regions... ";
        initRegions(mesh);
        std::cout << "done." << std::endl;
    }

    void
    BuiltInBuilder::buildMesh(getfem::mesh& mesh) const
    {

        size_type N = M_domain.dim;
        std::string meshType = M_domain.meshType; // element type for meshing, linear by defaults

        auto pgt = bgeot::geometric_trans_descriptor(meshType);

        if (true){
            auto pbDim = pgt->dim();
            assert(N==pbDim);
        }

        std::vector<size_type> nsubdiv(N);
        std::vector<scalar_type> lengths(N);
        nsubdiv[0] = M_domain.Nx;
        nsubdiv[1] = M_domain.Ny;
        nsubdiv[2] = M_domain.Nz;
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
    BuiltInBuilder::initRegions(getfem::mesh& mesh) const
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
                    // std::clog <<"["<<pt[0]<<","<<pt[1]<<","<<pt[2]<<"]\n";
                
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
                        // std::clog << "Adding an element in Fault convexList: "
                        // << "pt[0] = " << pt[0]
                        // << ", abs(pt[0]) = " << std::abs(pt[0]) << " < 1.e-6? "
                        // << std::boolalpha << (std::abs(pt[0]) < 1.e-6) << std::endl;
                    }
                    centeredConvexesList.add(i);
                    insertedCentered = true;
                }
            }
        }

        mesh.region(RegionType::BulkLeft).add(leftConvexesList);

        mesh.region(RegionType::BulkRight).add(rightConvexesList);

        for (dal::bv_visitor i(centeredConvexesList); !i.finished(); ++i){

            auto psc = mesh.structure_of_convex(i);

            // check the centroid (x < 0 || x > 0)      
            auto pts = mesh.points_of_convex(i);
            base_node centroid(pts[0].size());

            for (const auto& pt : pts)
                centroid += pt;
            centroid /= static_cast<double>(pts.size());

            // loop over the faces of a convex
            for (size_type f = 0; f < psc->nb_faces(); ++f){

                // check the normal
                auto pgt = mesh.trans_of_convex(i);
                auto face_pts = pgt->convex_ref()->points_of_face(f);

                base_small_vector n = mesh.mean_normal_of_face_of_convex(i, f);

                if (std::abs(n[0] - 1.0) < 1.e-6 && std::abs(n[1]) < 1.e-6 &&
                    std::abs(n[2]) < 1.e-6 && centroid[0] < -1.e-6){

                    // add the face to the fault's meshRegion
                    mesh.region(RegionType::Fault).add(i,f);
                }
            }
        }

        
        // Create the boundary regions
        getfem::mesh_region borderFaces;
        getfem::outer_faces_of_mesh(mesh, borderFaces);

        double eps = 1.e-2;

        auto make_box = [eps](const base_node& a, const base_node& b) {
            base_node pt_min(a.size()), pt_max(a.size());
            for (size_type i = 0; i < a.size(); ++i) {
                pt_min[i] = std::min(a[i], b[i]) - eps;
                pt_max[i] = std::max(a[i], b[i]) + eps;
            }
            return std::make_pair(pt_min, pt_max);
        };
    {
        auto [min1, max1] = make_box(
            base_node{-M_domain.Lx/2., 0., 0.},
            base_node{0., M_domain.Ly, 0.}
        );
        mesh.region(1) = getfem::select_faces_in_box(mesh, borderFaces, min1, max1);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{0., M_domain.Ly, 0.},
            base_node{-M_domain.Lx, M_domain.Ly, M_domain.Lz});
        mesh.region(2) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{-M_domain.Lx/2, M_domain.Ly, M_domain.Lz},
            base_node{0., 0., M_domain.Lz});
        mesh.region(3) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{0., 0., M_domain.Lz},
            base_node{-M_domain.Lx/2., 0., 0.});
        mesh.region(4) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }

    {       
        auto [pt_min, pt_max] = make_box(
            base_node{0., 0., 0.},
            base_node{M_domain.Lx/2., M_domain.Ly, 0.});
        mesh.region(5) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{M_domain.Lx/2., M_domain.Ly, 0.},
            base_node{0., M_domain.Ly, M_domain.Lz});
        mesh.region(6) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{0., M_domain.Ly, M_domain.Lz},
            base_node{M_domain.Lx/2., 0., M_domain.Lz});
        mesh.region(7) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {    
        auto [pt_min, pt_max] = make_box(
            base_node{M_domain.Lx/2., 0., M_domain.Lz},
            base_node{0., 0., 0.});
        mesh.region(8) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{-M_domain.Lx/2., 0., 0.},
            base_node{-M_domain.Lx/2., M_domain.Ly, M_domain.Lz});
        mesh.region(9) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }
    {
        auto [pt_min, pt_max] = make_box(
            base_node{M_domain.Lx/2., 0., 0.},
            base_node{M_domain.Lx/2., M_domain.Ly, M_domain.Lz});
        mesh.region(10) = getfem::select_faces_in_box(mesh, borderFaces, pt_min, pt_max);
    }

        // // alternatively:
        // // bds(9) = getfem::select_faces_of_normal(mesh, border_faces,
        // //    base_small_vector(-1.,0.,0.), 0.01);
        // // bds(10) = getfem::select_faces_of_normal(mesh, border_faces,
        // //    base_small_vector(1.,0.,0.), 0.01);

    }


    // GmshBuilder:: IMPLEMENTATION
    void
    GmshBuilder::buildMesh(getfem::mesh& mesh) const
    {
        using RegMap = std::map<std::string, size_type>;
        RegMap regmap;
        getfem::import_mesh_gmsh(M_meshFile, mesh, regmap);

        std::cout << regmap.size() << "\n";
        for (RegMap::iterator i = regmap.begin(); i != regmap.end(); i++) {
            std::cout << i->first << " " << i->second << "\n";
        }
    }

    void
    GmshBuilder::initRegions(getfem::mesh& mesh) const
    {
        getfem::mesh_region &fault_region = mesh.region(Fault);
        const getfem::mesh_region &bulk_region = mesh.region(BulkRight);

        // Loop through (convex, face) pairs
        std::clog << "Faces in region Fault (convex, face): ";
        for (getfem::mr_visitor it(fault_region); !it.finished(); ++it) {
            if (bulk_region.is_in(it.cv())){
                std::cout << "Removing (" << it.cv() << "," << it.f() << ")...";
                fault_region.sup(it.cv(),it.f());
            }
            std::cout << "done." << std::endl;
        }

    }

} // namespace gf