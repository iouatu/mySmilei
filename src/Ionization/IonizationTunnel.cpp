#include "IonizationTunnel.h"
#include "IonizationTables.h"

#include <cmath>

#include "Particles.h"
#include "Species.h"
using namespace std;



IonizationTunnel::IonizationTunnel( Params &params, Species *species ) : Ionization(params, species)
{
    DEBUG( "Creating the Tunnel Ionizaton class" );
    
    atomic_number_  = species -> atomic_number_;
    
    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    Potential.resize( atomic_number_ ); // Potential is a std::vector<double>, private member of IonizationTunnel class 
    Azimuthal_quantum_number.resize( atomic_number_ ); // Azimuthal_quantum_number is a std::vector<double>, private member of IonizationTunnel class
    
    for(int Zstar = 0; Zstar < (int)atomic_number_; Zstar++) {
        Potential               [Zstar] = IonizationTables::ionization_energy( atomic_number_, Zstar ) * eV_to_au; // returns the IonizEnergy of Element with atomic# = atomic_number_ and charge state Zstar, i.e. Element^{Zstar+}. Zstar=0 means neutral atom. conversion to AtomicUnits included at end.
        Azimuthal_quantum_number[Zstar] = IonizationTables::azimuthal_atomic_number( atomic_number_, Zstar );
    }
    
    for( unsigned int i=0; i<atomic_number_; i++ ) { // this for can be commented.
        DEBUG( "ioniz: i " << i << " potential: " << Potential[i] << " Az.q.num: " << Azimuthal_quantum_number[i] );
    } // this for can be commented.
    
    one_third = 1.0/3.0; // private member of IonizationTunnel class 

    alpha_tunnel.resize(atomic_number_); // alpha_tunel, beta_tunnel, gamma_tunnel
    beta_tunnel.resize(atomic_number_); //  are all std::vector<double>
    gamma_tunnel.resize(atomic_number_); // and are all private members of IonizationTunnel class.
    
    for( unsigned int Z=0 ; Z<atomic_number_ ; Z++ ) {
        DEBUG( "Z : " << Z );
        double cst      = ((double)Z + 1.0) * sqrt(2.0/Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z]  = pow(2, alpha_tunnel[Z]) * (8.*Azimuthal_quantum_number[Z] + 4.0) / (cst*tgamma(cst)) * Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0*Potential[Z], 1.5);
    }
    DEBUG( "Finished Creating the Tunnel Ionizaton class" );
}




void IonizationTunnel::operator()( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref )
{

    unsigned int Z, Zp1, newZ, k_times; // Z=0 is a neutral atom, Z=1 is a singly ionized atom. different than the Z above
    double TotalIonizPot, E, invE, factorJion, delta, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    vector<double> IonizRate_tunnel(atomic_number_), Dnom_tunnel(atomic_number_);
    LocalFields Jion;
    double factorJion_0 = au_to_mec2 * EC_to_au*EC_to_au * invdt; // Will be used to calc. ionization current.
    
    int nparts = Epart->size()/3; // Over 3 because there are 3 Dims. Epart is like [part1_Ex, part2_Ex, ... partn_Ex, part1_Ey, part2_Ey, ..., partn_Ey, part1_Ez, part2_Ez, ..., partn_Ez]
    double *Ex = &( (*Epart)[0*nparts] ); // Pointer to beginning of the vector which holds the (3 * n_parts) values of the n_parts E-field (Physics) 3D vectors.
    double *Ey = &( (*Epart)[1*nparts] ); // Pointer to where E_y values start inside the Epart vector.
    double *Ez = &( (*Epart)[2*nparts] ); // Pointer to where E_z values start inside the Epart vector. 
    
    for(unsigned int ipart=ipart_min ; ipart<ipart_max; ipart++) { // Loop on particles
    
        // Current charge state of the ion
        Z = (unsigned int) ( particles->charge(ipart) ); // Z=0 is a neutral atom. Z=1 is a singly ionized (positively charged) ion.
        
        // If ion already fully ionized then skip
        if( Z==atomic_number_ ) { continue; } // Hydrogen has atomic_number_=1. If Z=1 (singly ionized ion), the H ion is fully ionized. Makes sense.
        
        // Absolute value of the electric field normalized in atomic units
        E = EC_to_au * sqrt( pow( *( Ex + ipart - ipart_ref ), 2 ) // (Ex+ipart-ipart_ref) points to the location in the container Epart at which E-field_x for particle ipart sits.
                             +pow( *( Ey + ipart - ipart_ref ), 2 ) // Similarly for y. Dereferencing it via *() means we access and use the value sitting at that location.
                             +pow( *( Ez + ipart - ipart_ref ), 2 ) );
        if( E < 1e-10 ) {
            continue; // why is this?
        }
        
        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------
        
        invE = 1./E; // E is in AtomicUnits
        factorJion = factorJion_0 * invE*invE;
        delta      = gamma_tunnel[Z] * invE;
        ran_p = patch->rand_->uniform();
        IonizRate_tunnel[Z] = beta_tunnel[Z] * exp(-delta*one_third + alpha_tunnel[Z]*log(delta));
        
        // Total ionization potential (used to compute the ionization current)
        TotalIonizPot = 0.0;
        
        // k_times will give the number of ionization events
        k_times = 0;
        Zp1 = Z + 1;
        
        if( Zp1 == atomic_number_ ) { 
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if( ran_p < 1.0 - exp(-IonizRate_tunnel[Z]*dt) ) { // MC dictates to ionize it
                TotalIonizPot += Potential[Z];
                k_times        = 1; // ionize it
            } // else: nothing happens, no ionization of the last electron during this timestep. for loop moves to next particle, i.e. the one indexed by ipart+1.
        } 
        else {
            // multiple ionization can occur in one time-step
            //   partial & final ionization are decoupled (see Nuter Phys. Plasmas)
            // -------------------------------------------------------------------------
            
            // initialization
            Mult = 1.0;
            Dnom_tunnel[0] = 1.0;
            Pint_tunnel = exp( -IonizRate_tunnel[Z]*dt ); // cummulative prob. (Pint_tunnel is a double)
            
            // multiple ionization loop while Pint_tunnel < ran_p and still partial ionization
            while( (Pint_tunnel < ran_p) and (k_times < atomic_number_-Zp1) ) {
                newZ = Zp1 + k_times;
                delta = gamma_tunnel[newZ] * invE;
                IonizRate_tunnel[newZ] = beta_tunnel[newZ]
                                         * exp( -delta*one_third + alpha_tunnel[newZ]*log(delta));
                D_sum = 0.0;
                P_sum = 0.0;
                Mult  *= IonizRate_tunnel[ Z + k_times ];
                for(unsigned int i=0;  i < k_times+1; i++ ) {
                    Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z+i]);
                    D_sum += Dnom_tunnel[i];
                    P_sum += exp( -IonizRate_tunnel[Z+i]*dt ) * Dnom_tunnel[i];
                }
                Dnom_tunnel[k_times+1] -= D_sum;
                P_sum                   = P_sum + Dnom_tunnel[k_times+1] * exp(-IonizRate_tunnel[newZ]*dt);
                Pint_tunnel             = Pint_tunnel + P_sum*Mult;
                
                TotalIonizPot += Potential[ Z+k_times ];
                k_times++;
            } // END while
            
            // final ionization (of last electron)
            if( ((1.0-Pint_tunnel) > ran_p) && (k_times == atomic_number_-Zp1) ) {
                TotalIonizPot += Potential[atomic_number_-1];
                k_times++;
            } // END final ionization (of last electron)
        
        }// END Multiple ionization routine
        
        // Compute ionization current
        if (patch->EMfields->Jx_ != NULL) {  // For the moment ionization current is not accounted for in AM geometry
            factorJion *= TotalIonizPot;
            Jion.x = factorJion * *(Ex + ipart);
            Jion.y = factorJion * *(Ey + ipart);
            Jion.z = factorJion * *(Ez + ipart);
            
            Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, 
                                     patch->EMfields->Jz_, *particles, ipart, Jion);
        }  
        
        // Creation of the new electrons (variable weights are used)
        // -----------------------------
        if(k_times != 0) { // If ionization (single or multiple) has taken place.
            new_electrons.createParticle(); // new_electrons is an instance of class Particles.
            // 
            //new_electrons.initialize( new_electrons.size()+1, new_electrons.dimension() ); // ???
            int idNew = new_electrons.size() - 1; // new_electrons.size() =  number of electrons in new_electrons. -1 from idNew because the first elem of array has index 0.

            for(unsigned int i=0; i<new_electrons.dimension(); i++) { // new_electrons.dimension() is 1,2 or 3, depending on geometry of simulation.
                new_electrons.position(i, idNew) = particles->position(i, ipart);
            }
            for(unsigned int i=0; i<3; i++) { // Momentum is a 3D vector. (3-Velocity simulation)
                new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
            }

            new_electrons.weight(idNew) = double(k_times) * particles->weight(ipart); // 1e is generated, but its a special one. its weight is k_times higher than of simple 1e
            new_electrons.charge(idNew) = -1; // they are electrons with charge -e.
            
            // Increase the charge of the particle which ejected the electron(s).
            particles -> charge(ipart) += k_times; // however we add 1 e only to new_electrons object.

        } // END creation of electrons
        
    } // Loop on particles

} // void IonizationTunnel::operator() function's scope end
