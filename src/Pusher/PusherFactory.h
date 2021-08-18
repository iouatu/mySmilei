// --------------------------------------------------------------------------------------------------------------------
//
//! \file PusherFactory.h
//
//! \brief Class PusherFactory that manage the pusher association to the species.
//
// --------------------------------------------------------------------------------------------------------------------

#ifndef PUSHERFACTORY_H
#define PUSHERFACTORY_H

#include "Pusher.h"
#include "PusherBoris.h"
#include "PusherPonderomotiveBoris.h"
#include "PusherPonderomotivePositionBoris.h"
#include "PusherVay.h"
#include "PusherBorisNR.h"
#include "PusherRRLL.h"
#include "PusherHigueraCary.h"
#include "PusherPhoton.h"

#ifdef _VECTO
#include "PusherBorisV.h"
#include "PusherPonderomotiveBorisV.h"
#include "PusherPonderomotivePositionBorisV.h"
#endif

#include "Params.h"
#include "Species.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherFactory
//
//! \brief Class PusherFactory that manage the pusher association to the species.
//  --------------------------------------------------------------------------------------------------------------------
class PusherFactory
{
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate pusher for the species ispec
    //! \param ispec SpeciesId
    //! \param params Parameters
    //  --------------------------------------------------------------------------------------------------------------------
    static Pusher *create( Params &params, Species *species )
    {
        Pusher *Push = NULL;
        
        // Particle of matter
        if( species->mass_ > 0 ) {
            // assign the correct Pusher to Push
            if( species->pusher_name_ == "boris" ) {
                if( !species->vectorized_operators && !params.cell_sorting ) {
                    Push = new PusherBoris( params, species );
                }
#ifdef _VECTO
                else {
                    Push = new PusherBorisV( params, species );
                }
#endif
            } else if( species->pusher_name_ == "ponderomotive_boris" ) {
            
                int n_envlaser = params.Laser_Envelope_model;
                if( n_envlaser <1 ) {
                    ERROR( "No Laser Envelope present. The pusher ponderomotive_boris can be used only in presence of a Laser Envelope." );
                }
                
                if( !species->ponderomotive_dynamics ) {
                    ERROR( "if ponderomotive_boris pusher is chosen for a species, the flag ponderomotive_dynamics for that species must be set to true." );
                }
                
                if( !species->vectorized_operators && !params.cell_sorting ) {
                    Push = new PusherPonderomotiveBoris( params, species );
                }
#ifdef _VECTO
                else {
                    Push = new PusherPonderomotiveBorisV( params, species );
                }
#endif
            } else if( species->pusher_name_ == "borisnr" ) {
                Push = new PusherBorisNR( params, species );
            }
            /*else if ( species->pusher_name_ == "rrll" )
            {
                Push = new PusherRRLL( params, species );
            }*/
            else if( species->pusher_name_ == "vay" ) {
                Push = new PusherVay( params, species );
            } else if( species->pusher_name_ == "higueracary" ) {
                Push = new PusherHigueraCary( params, species );
            } else {
                ERROR( "For species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`" );
            }
        }
        // Photon
        else if( species->mass_ == 0 ) {
            if( species->pusher_name_ == "norm" ) {
                Push = new PusherPhoton( params, species );
            } else {
                ERROR( "For photon species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`" );
            }
        }
        
        if( species->ponderomotive_dynamics ) {
            if( species->pusher_name_ != "ponderomotive_boris" ) {
                ERROR( "For species " << species->name_ << " the flag ponderomotive_dynamics is True - the only pusher available to interact with the envelope is ponderomotive_boris" );
            }
        }
        return Push;
    }
    
    static Pusher *create_ponderomotive_position_updater( Params &params, Species *species )
    {
        Pusher *Push_ponderomotive_position = NULL;
        
        // Particle of matter
        if( species->mass_ > 0 ) {
            // assign the correct Pusher to Push_ponderomotive_position
            if( species->pusher_name_ == "ponderomotive_boris" ) {
                if( !species->vectorized_operators && !params.cell_sorting ) {
                    Push_ponderomotive_position = new PusherPonderomotivePositionBoris( params, species );
                }
#ifdef _VECTO
                else {
                    Push_ponderomotive_position = new PusherPonderomotivePositionBorisV( params, species );
                }
#endif
            }
            
            else {
                ERROR( "For species " << species->name_
                       << ": unknown pusher `"
                       << species->pusher_name_ << "`" );
            }
        } else {
            ERROR( "Ponderomotive pusher is not a valid choice for photons" );
        }
        return Push_ponderomotive_position;
    }
    
};

#endif
