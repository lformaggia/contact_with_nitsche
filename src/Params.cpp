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
        refined = command_line.search("-r");
        test = command_line.search("-t");
        init = command_line.search("-i");
        initRef = command_line.search("-ir");
        
        std::ifstream check(dataFileName);
        if(!check.is_open())
            throw std::runtime_error("Could not open the datafile!");
        check.close();
        
        domain.dim = datafile("domain/dim", 3);
        domain.Lx = datafile("domain/Lx", 1.0);
        domain.Ly = datafile("domain/Ly", 1.0);
        domain.Lz = datafile("domain/Lz", 1.0);
        domain.h = datafile("domain/h", 0.2);
        domain.angle = datafile("domain/angle", 0.);
        int ratio = static_cast<int>(domain.Ly/domain.Lx);
        domain.Nx = static_cast<int>(domain.Lx/domain.h);
        domain.Ny = static_cast<int>(domain.Ly/domain.h/ratio);
        domain.Nz = static_cast<int>(domain.Lz/domain.h/ratio);
        domain.meshType = datafile("domain/meshType", "GT_PK(3,1)");
        
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
        it.tol = datafile("it/tol", 1e-6);
        
        contact.method = datafile("contact/method", "nitsche");
        contact.theta = datafile("contact/theta", 0.0);
        contact.gammaN = datafile("contact/gammaN", 10.0)*physics.M_E0/domain.h; // gammaN = 10*E/h
        contact.gammaP = datafile("contact/gammaP", 10.0)*physics.M_E0/domain.h; // gammaP = 10*E/h
        contact.gammaL = datafile("contact/gammaL", 10.0)*physics.M_E0/domain.h; // gammaL = 10*E/h

        std::string order_u = datafile("domain/order_u", "1");
        std::string order_lm = datafile("domain/order_lm", "0");

        time.t0 = datafile("time/t0", 0.0);
        time.tend = datafile("time/tend", 1.0);
        time.dt = datafile("time/dt", 0.1);

        if (domain.meshType == "GT_PK(3,1)")
        {
            numerics.integration = "IM_TETRAHEDRON(5)";
            numerics.FEMTypeDisplacement = "FEM_PK(3,"+order_u+")";
            numerics.FEMTypeRhs = "FEM_PK(3,"+order_u+")";
            numerics.FEMTypeStress = "FEM_PK(3,"+order_u+")";
            numerics.FEMTypeLM = "FEM_PK(3,"+order_lm+")"; //"FEM_PK_DISCONTINUOUS(3,0)";
        }
        else if (domain.meshType == "GT_QK(3,1)")
        {
            numerics.integration = "IM_HEXAHEDRON(5)";
            numerics.FEMTypeDisplacement = "FEM_QK(3,"+order_u+")";
            numerics.FEMTypeRhs = "FEM_QK(3,"+order_u+")";
            numerics.FEMTypeStress = "FEM_QK(3,"+order_u+")";
            numerics.FEMTypeLM = "FEM_QK(3,"+order_lm+")";//"FEM_QK_DISCONTINUOUS(3,0)";
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
        os << "-- mu_friction = " << p.physics.M_mu_friction << "\n";
        os << "IT:\n";
        os << "--maxit = " << p.it.maxIt << "\n";
        os << "--tol = " << p.it.tol << "\n";
        os << "CONTACT: using "<< p.contact.method << " method\n";
        if (p.contact.method == "penalty") {
            os << "--gammaP = " << p.contact.gammaP << "\n";
        } else if (p.contact.method == "nitsche") {
            os << "--theta = " << p.contact.theta << "\n";
            os << "--gammaN = " << p.contact.gammaN << "\n";
        } else if (p.contact.method == "augLM"){
            os << "--gammaL = " << p.contact.gammaL << "\n";
        }
        else {
            throw std::runtime_error("Unsupported contact method: " + p.contact.method);
        }
        os << "TIME:\n";
        os << "--t0 = " << p.time.t0 << "\n";
        os << "--tend = " << p.time.tend << "\n";
        os << "--dt = " << p.time.dt << "\n";
        os << "NUMERICS:\n";
        os << "--integration = " << p.numerics.integration << "\n";
        os << "--FEMTypeDisplacement = " << p.numerics.FEMTypeDisplacement << "\n";
        os << "--FEMTypeStress = " << p.numerics.FEMTypeStress << "\n";
        os << "--FEMTypeRhs = " << p.numerics.FEMTypeRhs << "\n";
        if (p.contact.method == "augLM")
            os << "--FEMTypeLM = " << p.numerics.FEMTypeLM << "\n";
        os << "\nVERBOSE: " << (p.verbose ? "true" : "false") << "\n";
        os << "GMSH: " << (p.gmsh ? "true" : "false") << "\n";
        os << "============================================\n\n";
            
        return os;
    };

} // namespace gf