#include "Params.hpp"
#include "Utils.hpp"
#include "muParserXInterface.hpp"

namespace gf {

    Params::Params(int argc, char* argv[])
    : datafile("data.pot")
    {
        GetPot command_line(argc, argv);
        const std::string dataFileName = command_line.follow("data.pot", 2, "-f",
                "--file");
        
        verbose = command_line.search("-v");
        gmsh = command_line.search("-m");
        
        std::ifstream check(dataFileName);
        if(!check.is_open())
            throw std::runtime_error("Could not open the datafile!");
        check.close();
        
        // GetPot datafile {filename.c_str()};
        
        domain.dim = datafile("domain/dim", 3);
        domain.Lx = datafile("domain/Lx", 1.0);
        domain.Ly = datafile("domain/Ly", 1.0);
        domain.Lz = datafile("domain/Lz", 1.0);
        domain.h = datafile("domain/h", 0.2);
        domain.angle = datafile("domain/angle", 0.);
        domain.Nx = static_cast<int>(domain.Lx/domain.h);
        domain.Ny = static_cast<int>(domain.Ly/domain.h);
        domain.Nz = static_cast<int>(domain.Lz/domain.h);
        domain.meshType = datafile("domain/meshType", "GT_QK(3,1)");
        
        physics.M_E0 = datafile("physics/E", 0.0);
        physics.M_nu = datafile("physics/nu", 0.0);
        physics.M_mu_friction = datafile("physics/mu_friction", 0.5);
        physics.M_lambda = (physics.M_E0 * physics.M_nu) /
                            ((1 + physics.M_nu) * (1 - 2 * physics.M_nu));
        physics.M_mu = physics.M_E0 / (2 * (1 + physics.M_nu));
        // Load gravity vector
        std::string gravity_str = datafile("physics/bulkLoad", "[0., 0., 0.]"); // Note: expecting [ ... ]
        std::vector<std::string> gVecStr = splitString(gravity_str);
        for (size_type i{}; i < 3; ++i)
            physics.M_gravity.push_back(std::stod(gVecStr[i]));

        if (physics.M_gravity.size() != 3)
            throw std::runtime_error("Expected 3 components in physics/bulkLoad");

        it.maxIt = datafile("it/maxit", 30);
        it.rtol = datafile("it/tol", 1e-6);
        it.atol = datafile("it/atol", 0.); // atol: not in file, default to 0
        
        nitsche.theta = datafile("it/theta", 0.0);
        nitsche.gamma0 = datafile("it/gamma", 10.0);

        time.t0 = datafile("time/t0", 0.0);
        time.tend = datafile("time/tend", 1.0);
        time.dt = datafile("time/dt", 0.1);

        if (domain.meshType == "GT_PK(3,1)")
        {
            numerics.integration = "IM_TETRAHEDRON(5)";
            numerics.FEMTypeDisplacement = "FEM_PK(3,1)";
            numerics.FEMTypeRhs = "FEM_PK(3,1)";
            numerics.FEMTypeStress = "FEM_PK(3,1)";
        }
        else if (domain.meshType == "GT_QK(3,1)")
        {
            numerics.integration = "IM_HEXAHEDRON(5)";
            numerics.FEMTypeDisplacement = "FEM_QK(3,1)";
            numerics.FEMTypeRhs = "FEM_QK(3,1)";
            numerics.FEMTypeStress = "FEM_QK(3,1)";
            // std::cout << "Inside Params constructor: numerics.FEMTypeRhs = " << numerics.FEMTypeRhs << std::endl;
        }
        else 
            throw std::runtime_error("Select either GT_PK(3,1) or GT_QK(3,1)");

        if (verbose) std::cout << *this << std::endl;
    }

    std::ostream& operator<<(std::ostream& os, const Params& p){
        os << "================ PARAMETERS ================\n";
        os << "DOMAIN:\n";
        os << "-- [Lx, Ly, Lz] = ["
           << p.domain.Lx << ", "
           << p.domain.Ly << ", "
           << p.domain.Lz << "]\n";
        os << "-- [Nx, Ny, Nz] = ["
           << p.domain.Nx << ", "
           << p.domain.Ny << ", "
           << p.domain.Nz << "]\n";
        os << "PHYSICS:\n";
        os << "-- E0 = " << p.physics.M_E0 << "\n";
        os << "-- nu = " << p.physics.M_nu << "\n";
        os << "-- g = [" 
           << p.physics.M_gravity[0] << ", "
           << p.physics.M_gravity[1] << ", "
           << p.physics.M_gravity[2] << "]\n";
        os << "IT:\n";
        os << "--maxit = " << p.it.maxIt << "\n";
        os << "--rtol = " << p.it.rtol << "\n";
        os << "--atol = " << p.it.atol << "\n";
        /** \todo: print other parameters */
    
        return os;
    };

} // namespace gf