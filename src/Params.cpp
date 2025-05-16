#include "Params.hpp"
#include "Utils.hpp"
#include "muParserXInterface.hpp"

namespace gf {

    Params::Params(const std::string& filename, const std::string& meshfile, bool ver)
        : datafile{
            filename.c_str()
        },
        domain{
            datafile("domain/dim", 3),
            datafile("domain/Lx", 1.0),
            datafile("domain/Ly", 1.0),
            datafile("domain/Lz", 1.0),
            datafile("domain/Nx", 1),
            datafile("domain/Ny", 1),
            datafile("domain/Nz", 1),
            datafile("domain/meshType", "GT_QK(3,1)")
        },
        physics{
            datafile("physics/E", 0.0),
            datafile("physics/nu", 0.0),
            0.0, // M_lambda (computed below)
            0.0, // M_mu (computed below)
            {}   // M_gravity (set below)
        },
        it{
            datafile("it/ContactMaxIt", 100),
            datafile("it/ContactToll", 1e-6),
            0.0 // atol: not in file, default to 0
        },
        nitsche{
            datafile("it/theta", 0.0),
            datafile("it/gamma", 10.0)
        },
        time{
            datafile("time/t0", 0.0),
            datafile("time/tend", 1.0),
            datafile("time/dt", 0.1)
        },
        numerics{
            datafile("numerics/integration","IM_HEXAHEDRON(5)")
        },
        meshFile{meshfile},
        verbose{ver}
   
    {
        // Compute dependent physics values
        physics.M_lambda = (physics.M_E0 * physics.M_nu) /
                            ((1 + physics.M_nu) * (1 - 2 * physics.M_nu));
        physics.M_mu = physics.M_E0 / (2 * (1 + physics.M_nu));
    
        // Load gravity vector
        std::string gravity_str = datafile("physics/bulkLoad", "[0.0 0.0 0.]"); // Note: expecting [ ... ]
        std::vector<std::string> gVecStr = splitString(gravity_str);
        for (size_type i{}; i < 3; ++i)
            physics.M_gravity.push_back(std::stod(gVecStr[i]));

        if (physics.M_gravity.size() != 3)
            throw std::runtime_error("Expected 3 components in physics/bulkLoad");

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