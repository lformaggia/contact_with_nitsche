#include "DomainView.hpp"
#include "Domain.hpp"

namespace gf{
    
    FaultView::FaultView(getfem::mesh& mesh, const dal::bit_vector& centeredConvexesList)
    : DomainView(mesh)
    {
        // Add the FACES of the convexes s.t. they belong to the interface
        // M_normal = /*...*/;
        // loop  over the convexes
        for (dal::bv_visitor i(centeredConvexesList); !i.finished(); ++i){

            auto psc = M_mesh.structure_of_convex(i);

            // check the centroid (x < 0 || x > 0)      
            auto pts = M_mesh.points_of_convex(i);
            base_node centroid(pts[0].size());

            for (const auto& pt : pts)
                centroid += pt;
            centroid /= static_cast<double>(pts.size());

            // loop over the faces of a convex
            for (size_type f = 0; f < psc->nb_faces(); ++f){

                // check the normal
                auto pgt = M_mesh.trans_of_convex(i);
                auto face_pts = pgt->convex_ref()->points_of_face(f);

                base_small_vector n = M_mesh.mean_normal_of_face_of_convex(i, f);
                if (std::abs(n[0] - 1.0) < 1e-6 && std::abs(n[1]) < 1e-6 && std::abs(n[2]) < 1e-6){
                    // add the face to the fault's meshRegion
                    
                    M_region.add(f);

                    // take the neighbor element through face f

                    size_type neigh = M_mesh.neighbor_of_convex(i,f);

                    size_type right = (centroid[0] < -1e-6) ? i : neigh;
                    size_type left = (right == i) ? neigh : i;

                    // add the pair (left,right) to a map
                    convexesSharingFace[f] = std::make_pair(left,right);

                }
            }
        }
    }

}
