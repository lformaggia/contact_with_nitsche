#include "MeshRegion.hpp"

bool DEBUGMR = false;

namespace gf {
    
    // Fault::Fault(const getfem::mesh& mesh, const dal::bit_vector& centeredConvexesList)
    // : MeshRegion(mesh)
    // {
    //     // Add the FACES of the convexes s.t. they belong to the interface
    //     // M_normal = /*...*/;
    //     // loop  over the convexes
    //     for (dal::bv_visitor i(centeredConvexesList); !i.finished(); ++i){

    //         auto psc = M_mesh.structure_of_convex(i);

    //         // check the centroid (x < 0 || x > 0)      
    //         auto pts = M_mesh.points_of_convex(i);
    //         base_node centroid(pts[0].size());

    //         for (const auto& pt : pts)
    //             centroid += pt;
    //         centroid /= static_cast<double>(pts.size());

    //         if (DEBUGMR)
    //             std::clog << "Element "<< i
    //                 << ": x_c = ["<<centroid[0]
    //                 <<","<<centroid[1]<<","<<centroid[2]<<"]"
    //                 << std::endl;

    //         // loop over the faces of a convex
    //         for (size_type f = 0; f < psc->nb_faces(); ++f){

    //             // check the normal
    //             auto pgt = M_mesh.trans_of_convex(i);
    //             auto face_pts = pgt->convex_ref()->points_of_face(f);

    //             base_small_vector n = M_mesh.mean_normal_of_face_of_convex(i, f);

    //             if (DEBUGMR)
    //                 std::clog << "-- Face "<<f<<": n = ["
    //                     <<n[0]<<","<<n[1]<<","<<n[2]<<"]";

    //             if (std::abs(n[0] - 1.0) < 1.e-6 && std::abs(n[1]) < 1.e-6 &&
    //                 std::abs(n[2]) < 1.e-6 && centroid[0] < -1.e-6){
    //                 // add the face to the fault's meshRegion
    //                 if (DEBUGMR)
    //                     std::clog << " --> Adding to the fault region!\n";
                    
    //                 /** @remark this is actually unnecessary, since the addition
    //                  * of a face in a mesh region adds the index of the convex with
    //                  * the corresponding local face
    //                  * Since I decided to add the faces from the top (where i want to compute inegrals),
    //                  * the top element with its local face is added to the region.
    //                  * Thus, when looping over the Fault region, to get the neighboring element just use
    //                  * the left element idx with the local face, i.e. i.cv(), i.f()
    //                  * @todo: maybe the convexSharingFacesMap is unnecessary, remove if not needed
    //                  */
    //                 M_region.add(i,f);
    //                 getfem::mr_visitor j(M_region);

    //                 for (size_type jj{}; jj < M_region.nb_convex() - 1;++jj)
    //                     ++j;

    //                 // take the neighbor element through face f
    //                 size_type neigh = M_mesh.neighbor_of_convex(i,f);

    //                 size_type right = (centroid[0] < -1e-6) ? i : neigh;
    //                 size_type left = (right == i) ? neigh : i;

    //                 // add the pair (left,right) to a map

    //                 if (DEBUGMR){
    //                     std::clog << "Making pair for face "<<f<<std::endl;
    //                 }
    //                 convexesSharingFace[j.cv()] = std::make_pair(left,right);
    //                 ++j;
    //             }
    //             if (DEBUGMR)
    //                 std::clog << std::endl;


    //         }
    //     }
    //     if (DEBUGMR){
    //         std::clog << "Elements in the Fault region: ";
    //         for (getfem::mr_visitor i(M_region); !i.finished(); ++i){
    //             std::clog << i.cv() << " - face ";
    //             std::clog << i.f() << "\n";
    //         }
    //         std::clog << std::endl;

    //         for (const auto& [f,cv]: convexesSharingFace)
    //             std::clog << "Face "<<f
    //                 <<" is between elements ("
    //                 <<cv.first<<","<<cv.second<<")\n";
    //         std::clog << std::endl;
    //     }
    // }

}
