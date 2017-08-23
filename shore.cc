/******************************************************************************
 $Id: flux_stationary.cc,v 1.27 2005/04/13 18:33:12 duran Exp $
 ******************************************************************************/

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <iostream>

#include "physics_const.h"
#include "globals.h"

#include "shore.h"


//*****************************************************************************
//  class flux3d_stationary

/*!  Constructor: reads relevant parameters and creates auxiliary arrays for
 the velocity and the density of entrained sand rho and its gradient.  */

shore3d::shore3d( const dunepar& parameters ): dunedata(parameters)
{
    
    m_x_periodic= duneglobals::periodic_x();
    m_y_periodic= duneglobals::periodic_y();
    
    //!! BEACH
    m_shore_HMWL = duneglobals::HMWL(); //parameters.getdefault("shore.HMWL", 0.0);
    m_watertable = duneglobals::MSL(); //parameters.getdefault("shore.sealevel", 0.0);
    m_slope = duneglobals::slope();//parameters.getdefault("beach.angle", 0.0);
        
    m_sealevelrise = parameters.getdefault("shore.sealevelrise", 0.0); // m/yr
    m_sealevelrise /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec
    
//    m_shoreface_lenght = parameters.getdefault("shore.facelenght", 1000.0);
    m_grad_alongshore = parameters.getdefault("shore.alongshore_grad", 0.0); // m/yr
    m_grad_alongshore *= m_slope;
    m_grad_alongshore /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec

    m_shoreline = 0;

    // 
    
}


/*!  The destructor currently does nothing.  */
shore3d::~shore3d()
{
}

void shore3d::shorelinecal(const TFktScal& h)
{
    // CALCULATION OF THE MIN SHORELINE POSITION
    double shorelinepos = duneglobals::nx();
    int xaux = 0;
    for(int y = 0; y< duneglobals::ny(); ++y ){
        for(int x = 0; x< duneglobals::nx(); ++x ){
            if (h(x, y) >= m_shore_HMWL){
                xaux = x;
                break;
            }
        }
        // take smaller position
        shorelinepos = (shorelinepos > xaux ? xaux : shorelinepos);
    }   
    m_shoreline = shorelinepos; // store shoreline position
}

// correct shoreline geometry
void shore3d::restoreshoreface(TFktScal& h)
{
    int x0 = m_shoreline-5;
    for( int y= 0; y< duneglobals::ny(); ++y ){
        for( int x= x0; x < m_shoreline+1; ++x ){
            h(x, y) = h(x0, y) + m_slope*duneglobals::dx()*(x-x0);
        }
    }
}

int shore3d::shorefacemotion(TFktScal& h, double timestep)
{
   
    int mean_shoreshift = 0;
    double shoreshift_rate = m_sealevelrise + m_grad_alongshore;
    
    // CONSTANT SLOPE!
    int xend = (m_sealevelrise > 0 ? duneglobals::nx() : m_shoreline);
    for( int y= 0; y< duneglobals::ny(); ++y ){
        for( int x= 0; x < xend; ++x ){
            
            h(x, y) += - shoreshift_rate * timestep; // * (h(x, y) < m_shore_HMWL ? 1 : 1);
            
        }
        if (shoreshift_rate > 0) {
            // EROSION
            mean_shoreshift = (h(1,0) < 0 ? 1 : 0);
            
            // correct shoreline geometry
            restoreshoreface(h);
            
        } else {
            // ACCRETION
            mean_shoreshift = (h(0,0) >= m_slope * duneglobals::dx() ? -1 : 0);
        }
    }
    // CORRECT SHORELINE
    
    double shift = h(0,0)/m_slope;
        
    cout << h(0,0) << " " << m_shoreline << " " << m_slope << ' ' << mean_shoreshift << " # SHORE 1" << endl;
    
    return mean_shoreshift;
}


/*!  Saves the arrays  */

void shore3d::save_arrays()
{
}

