/******************************************************************************
 $Id: flowbeach.h,v 1.5 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#ifndef STORM_H
#define STORM_H

#include "globals.h"
#include "func.h"
#include "PTG_FileNames.h"
#include "vec.h"

class dunepar;
class avalanche;

/*! Continuum aproach of sand relaxation by coasts (including coasts over solid substrates) */

class storm : public dunedata
{

public:
    // construction
    storm(const dunepar& p);
    virtual ~storm() {}
        
    void CalcGradUp(TFktScal &h);
    void Step(TFktScal &h, TFktScal &overwash);

    virtual double impact(double shoreline, TFktScal &h, TFktScal &h_nonerod, TFktScal &overwash );
    virtual void stop( double time, double timestep, bool &calc_storm);

    virtual void calc( TFktScal &h, TFktScal &overwash );
  
    virtual void save_arrays();
        
private:
 
    /*!  Avalanche relaxation object.  */
    avalanche *m_avalanche;

    TFktScal m_sflux;
    TFktScal m_hst;
    
    int m_storm_iter;
    double m_Smax;
    double m_Sdt;
    double m_Q, m_scalefactor, m_shore_HMWL, m_watertable, m_surge, m_Tsurge;
    double m_slope;
    
    int m_shoreline;
    
    double surge[3000];
    int stormindex;

};


#endif

