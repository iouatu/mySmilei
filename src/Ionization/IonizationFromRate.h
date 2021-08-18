#ifndef IONIZATIONFROMRATE_H
#define IONIZATIONFROMRATE_H

#include <cmath>

#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;

//! calculate the particle FromRate ionization
class IonizationFromRate : public Ionization
{
public:
    IonizationFromRate( Params &params, Species *species ); //! CTOR for IonizationFromRate: with no input argument. Will be overriden in .cpp
    
    //! apply the FromRate Ionization model to the species
    void operator()( Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ) override;
    //! () above overrides the virtual operator() function from Ionization class.

private:
    int itime;
    unsigned int maximum_charge_state_;
    PyObject *ionization_rate_;
};

#endif
