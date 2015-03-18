/******************************************************************************
 $Id: flux_stationary.h,v 1.14 2005/04/13 17:47:53 duran Exp $
 ******************************************************************************/

#ifndef SHORE_H
#define SHORE_H

#include "vec.h"
#include "func.h"
#include "globals.h"

class dunepar;

/*!  \brief Class defining the main shore processes.  */

class shore3d : public dunedata
{
public:
    shore3d( const dunepar& parameters );
    virtual ~shore3d();

    /*!  Returns m_shoreline.  */
    virtual double shorelinepos() { return m_shoreline; }

    virtual void shorelinecal(const TFktScal& h);
    
    virtual void save_arrays();
    
protected:
    
    bool m_x_periodic, m_y_periodic;
    
    double m_watertable, m_shore_HMWL, m_slope, m_shoreline;


};


#endif // SHORE_H
