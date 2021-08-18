#include "IonizationTables.h"

// Gets the ionization energy of a given ion
double IonizationTables::ionization_energy(int atomic_number, int Zstar)
{
    return ionizationEnergy[( atomic_number*(atomic_number-1) )/2 + Zstar]; // returns the ionizEnergy of element with atomic number = atomic_number_ and charge state Zstar (i.e. Zstar=0 is a neutral atom He^{0}, Zstar=1 is singly ionized atom, i.e. Na^{+1}, etc.)
}


// Gets the azimuthal atomic number of a given ion
int IonizationTables::azimuthal_atomic_number(int atomic_number, int Zstar)
{
    return (int)azimuthalQuantumNumber[( atomic_number*( atomic_number-1 ) )/2 + Zstar]; // why this index? why this index? returns the ionizEnergy of element with atomic number = atomic_number_ and charge state Zstar (i.e. Zstar=0 is a neutral atom He^{0}, Zstar=1 is singly ionized atom, i.e. Na^{+1}, etc.)
}


// Gets the k-th binding energy in any neutral or ionized atom with atomic number Z and charge Zstar
// We use the formula by Carlson et al., At. Data Nucl. Data Tables 2, 63 (1970)
double IonizationTables::binding_energy(int atomic_number, int Zstar, int k)
{
    int offset = ( atomic_number * (atomic_number-1) ) / 2;
    return (( ionizationEnergy[offset + Zstar]
             - bindingEnergy   [offset + atomic_number - Zstar-1]
             + bindingEnergy   [offset + k                    ]
           ) / 510998.9); // converted to mc^2
}

