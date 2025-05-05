#include "Params.hpp"

namespace gf {

    Params::Params(const GetPot& datafile)
        : domain{
            datafile("domain/Lx", 1.0),
            datafile("domain/Ly", 1.0),
            datafile("domain/Lz", 1.0),
            datafile("domain/Nx", 1),
            datafile("domain/Ny", 1),
            datafile("domain/Nz", 1)
        },
            phyisics{
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
            }
    {
        // Compute dependent physics values
        phyisics.M_lambda = (phyisics.M_E * phyisics.M_nu) /
                            ((1 + phyisics.M_nu) * (1 - 2 * phyisics.M_nu));
        phyisics.M_mu = phyisics.M_E / (2 * (1 + phyisics.M_nu));
    
        // Load gravity vector
        phyisics.M_gravity = datafile("physics/bulkLoad", base_small_vector{0.0, 0.0, 0.0});
        if (g.size() != 3)
            throw std::runtime_error("Expected 3 components in physics/bulkLoad");

    }
    
} // namespace gf