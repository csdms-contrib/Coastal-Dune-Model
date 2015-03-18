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
        
    m_shoreline = 0;
    
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

/*!  Saves the arrays  */

void shore3d::save_arrays()
{
}

