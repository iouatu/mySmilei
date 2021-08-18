#ifndef IONIZATIONTUNNELBSI_H
#define IONIZATIONTUNNELBSI_H

#include <cmath>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;


class IonizationTunnelBSI : public Ionization {
    public:
        IonizationTunnelBSI(Params &params, Species *species); // C-TOR
        void operator()(Particles*, unsigned int, unsigned int, std::vector<double>*, Patch*, Projector*, int ipart_ref = 0) override; // will be ovverriden in the .cpp file assoc with this .h

    private:
        unsigned int atomic_number_; // the atomic # of the species. obtained from an input in the .py namelist the user provides. 
        std::vector<double> Potential, Azimuthal_quantum_number; // Potential contains IonizationPotentials in AtomicUnits for each charge state Z=0, 1,... atomic_number_. Z=0 is neutral atom.
        double one_third;
        std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel; // for the Tunnel ADK-PPT rate formula

};

#endif