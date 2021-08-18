#ifndef IONIZATIONTUNNEL_H
#define IONIZATIONTUNNEL_H

#include <cmath>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;


class IonizationTunnel : public Ionization
{

public:
    //! Constructor for IonizationTunnel: with no input argument
    IonizationTunnel( Params &params, Species *species );
    
    //! apply the Tunnel Ionization model to the species (with ionization current)
    void operator()( Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ) override;
    
private:
    unsigned int atomic_number_;
    std::vector<double> Potential;
    std::vector<double> Azimuthal_quantum_number;
    double one_third;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel;
};



#endif
