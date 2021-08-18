#ifndef IONIZATION_H
#define IONIZATION_H

#include <map>

#include "Tools.h"
#include "Params.h"
#include "Patch.h"
#include "Field.h"
#include "Particles.h"
#include "Projector.h"


//! Class Ionization: generic class allowing to define Ionization physics
class Ionization
{

public:
    Ionization( Params &params, Species *species ); // C-TOR of Ionization class.
    virtual ~Ionization(); // ~ signifies the D-TOR of Ionization class.
    
    //! Overloading of () operator. Defined as virtual to be overriden in child classes (i.e. IonizationFromRate, etc.)
    virtual void operator()( Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ) {};
    //! method for envelope ionization
    virtual void envelopeIonization( Particles *, unsigned int, unsigned int, std::vector<double> *, std::vector<double> *, std::vector<double> *, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ){};

    Particles new_electrons; // Particles is a class. Its C-TOR calls: short_prop.resize(0).
    
protected:

    double eV_to_au;
    double au_to_mec2;
    double EC_to_au; // helps normalize E-fields to atomic electric field unit (i.e. transforms E-fields from SI to atomic units) 
    double au_to_w0; 
    
    double reference_angular_frequency_SI;
    double dt;
    double invdt;
    unsigned int nDim_field;
    unsigned int nDim_particle;
    double ionized_species_invmass;
    
private:


};

#endif
