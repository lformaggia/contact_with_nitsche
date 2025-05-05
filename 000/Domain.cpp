#include "Domain.hpp"

namespace gf
{
    // Domain:: IMPLEMENTATION

    Domain::Domain(const GetPot& datafile) :
    M_mu(datafile("physics/mu", 1.0)),
    M_lambda(datafile("physics/lambda", 1.0))
    {
        // check if any mesh file is given
        std::string externalMesh = datafile("domain/meshFile", "");
        if (!externalMesh.empty())
            M_meshBuilder = std::make_unique<GmshBuilder>(externalMesh);
        else
            M_meshBuilder = std::make_unique<BuiltInBuilder>(datafile);

        // build the mesh
        M_meshBuilder->buildMesh(*this);

        // initialize the regions
        M_meshBuilder->initRegions(*this);
    }


    void Domain::exportMesh(const std::string &filename) const {
        getfem::vtk_export vtkmesh(filename);
        vtkmesh.exporting(M_mesh);
    }


    // BuiltInBuilder:: IMPLEMENTATION    

    void BuiltInBuilder::buildMesh(Domain& domain) const
    {
        size_type N = M_datafile("domain/problemDimension", 3);
        size_type Lx = M_datafile("domain/Lx", 2); // lenghts
        size_type Ly = M_datafile("domain/Ly", 4);
        size_type Lz = M_datafile("domain/Lz", 4);
        size_type Nx = M_datafile("domain/Nx", 4); // number of subdivisions
        size_type Ny = M_datafile("domain/Ny", 16);
        size_type Nz = M_datafile("domain/Nz", 16);
        std::string meshType = M_datafile("domain/meshType",
            (N < 3) ? "GT_PK(3,1)" : "GT_PK(3,1)"); // element type for meshing, linear by defaults


        auto pgt = bgeot::geometric_trans_descriptor(meshType);
        if (true){
            auto pbDim = pgt->dim();
            assert(N==pbDim);
        }

        std::vector<size_type> nsubdiv(N);
        std::vector<size_type> lengths(N);
        nsubdiv[0] = Nx;
        nsubdiv[1] = Ny;
        lengths[0] = Lx;
        lengths[1] = Ly;

        if (N==3)
        {
            nsubdiv[2] = Nz;
            lengths[2] = Lz;
        }
    
        // create the mesh
        getfem::regular_unit_mesh(domain.M_mesh, nsubdiv, pgt);

        // transform the mesh (scale the unit mesh to [Lx, Ly(, Lz)])
        //!\todo: maybe I also need to translate it, I want the fault to be positioned in x=0
        bgeot::base_matrix M(N,N);
        for (size_type i = 0; i < N; ++i)
            M(i,i) = lengths[i];

        domain.M_mesh.transformation(M);
        
        //!\todo: need to be sure that i can build the fracture according to the normal ...

    }

    void BuiltInBuilder::initRegions(Domain &domain) const {
        
        dal::bit_vector leftConvexesList;
        dal::bit_vector rightConvexesList;
        dal::bit_vector centeredConvexesList;
        
        // loop over all convexes in mesh
        for (dal::bv_visitor i(domain.M_mesh.convex_index()); !i.finished(); ++i){
            
            bool insertedLeft = false, insertedRight = false, insertedCentered = false;

            // get the coords of the points of convex i
            auto local_points = domain.M_mesh.points_of_convex(i);

            // loop over the points of a convex
            for (const auto &pt: local_points){
                if (pt[0] < -1e-6 && !insertedLeft){ // if one point has x < 0 --> left bulk region
                    leftConvexesList.add(i);
                    insertedLeft = true;
                }
                else if (pt[0] > 1.e-6 && !insertedRight){ // "   "   "   x > 0 --> right   "   "
                    rightConvexesList.add(i);
                    insertedRight = true;
                }
                else if (!insertedCentered)
                    centeredConvexesList.add(i);
                    insertedCentered = true;
            }
        }
        domain.M_regions["BulkLeft"] = std::make_unique<BulkView>(domain.M_mesh, leftConvexesList, SideType::LEFT);
        domain.M_regions["BulkRight"] = std::make_unique<BulkView>(domain.M_mesh, rightConvexesList, SideType::RIGHT);
        domain.M_regions["Fault"] = std::make_unique<FaultView>(domain.M_mesh,centeredConvexesList);
    }


    // GmshBuilder:: IMPLEMENTATION 

    void GmshBuilder::buildMesh(Domain &domain) const
    {
        domain.M_mesh.read_from_file(M_meshFile);
    }

    void GmshBuilder::initRegions(Domain &domain) const {
        std::map<std::string, size_type> regionMap;
        getfem::import_mesh_gmsh(M_meshFile, domain.M_mesh, regionMap);

        /* ... */
    
        auto itLeft = regionMap.find("BulkLeft");
        auto itRight = regionMap.find("BulkRight");
        auto itFault = regionMap.find("Fault");
    
        if (itLeft != regionMap.end()) {
            auto view = std::make_unique<BulkView>(domain.M_mesh, itLeft->second, SideType::LEFT);
            domain.M_regions["BulkLeft"] = std::move(view);
        }
    
        if (itRight != regionMap.end()) {
            auto view = std::make_unique<BulkView>(domain.M_mesh, itRight->second, SideType::RIGHT);
            domain.M_regions["BulkRight"] = std::move(view);
        }
    
        if (itFault != regionMap.end()) {
            auto view = std::make_unique<FaultView>(domain.M_mesh, itFault->second);
            domain.M_regions["Fault"] = std::move(view);
        }
    }
    
}