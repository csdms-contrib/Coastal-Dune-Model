/******************************************************************************
 $Id: vegetation.cc,v 1.9 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#include <math.h>

#include "globals.h"
#include "func.h"
#include "initsurf.h"
#include "vegetation.h"

//*****************************************************************************
//  class vegetation

vegetation::vegetation(const dunepar& p) : dunedata(p)
{
    m_veget.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );
    
    /*Vegetation parameters*/
    m_Lveg = p.getdefault("veget.xmin", 0.0);
    m_xmin0 = m_Lveg/duneglobals::dx();
    m_xmin = m_xmin0;

    m_Tveg = p.getdefault("veget.Tveg", 1.0) * duneglobals::secday() * duneglobals::timefrac();       // (sec)
    cout << "grass constructor: Tveg = " << m_Tveg << endl;
        
    // erosion/acc
    m_sens = p.getdefault("veget.sensitivity", 1.0);       // species sensitivity for erosion/acc. rate
    m_Hveg = p.getdefault("veget.Hveg", 1.0);       // maximum plant height (~ 1m)
    m_sens /= m_Hveg;
    
    // time conversion
    m_wind_factor =  duneglobals::timefrac();
        
    // extra
    m_rho_max = p.getdefault("veget.rho.max", 1.0);    
    
}

/*! Initialize vegetation */
void vegetation::init(const dunepar& par)
{
    arrayinit *init_veget;
    
    if( par.exists("veget.Init-Surf") )
        init_veget= arrayinit::create(par, "veget.");
    else
        init_veget= new CInitSurfPlain(0.0);
    
    init_veget->init_2d_vec( m_veget );
    
    delete init_veget;
    
}

/*! return cover fraction*/
void vegetation::getcover(TFktScal& rho_veget)
{
    // NORMALIZATION
    for( int x= 0; x< duneglobals::nx(); ++x )
        for( int y= 0; y< duneglobals::ny(); ++y ){
            // added density:
            rho_veget(x,y) = m_veget(x,y)[0] + 0*m_veget(x,y)[1];
            if (rho_veget(x,y) > m_rho_max) {
                rho_veget(x,y) = m_rho_max;
            }
        }
    
}

/*!  Computes the evolution of the vegetation.  */
int vegetation::evol(TFktScal& rho_veget, const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dh_dt)
{
    int veget_X0 = 0;
    
    veget_X0 = evolspec(time, timestep, shoreline, h, dh_dt, 0);
    
    getcover(rho_veget);
    
    return veget_X0;
}

/*!  Computes the evolution of each species.  */
int vegetation::evolspec(const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dh_dt, int species)
{
    // calculate gradient
    TFktVec grad_h;
    grad_h.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h.GradMid(h);

    // General variables
    double growthrate, proprate, reprod_rate, mortal_rate;
    int vegetpoints = 0;
    
    double m_shore_HMWL = duneglobals::HMWL();
    
    
    // Species growth parameters:
    // GENERIC (1)
    double V_gen = 1./m_Tveg;
        
    for(int y = 0; y< duneglobals::ny(); ++y ){
        for(int x = 0; x< duneglobals::nx()-1; ++x ){
            
            // AUXILIAR
            double dhdt = dh_dt(x,y); // erosion rate
            double erosion = (dhdt < 0 ? 1 : 0);
                        
            // define an alternative cover density encoding competition
            double rho_competition = m_veget(x,y)[0];
            rho_competition = (rho_competition > m_rho_max ? m_rho_max : rho_competition);
            // species limiting factors
            double shorefactor = (x < shoreline + m_xmin ? 0 : 1);  // 0: maximum effect; 1: no effect
          
            // GENERIC GROWTH
            double V = (1 - rho_competition) * V_gen * shorefactor;
            
            reprod_rate = V;
            mortal_rate = erosion * fabs(dhdt) * m_sens;
            
            growthrate = reprod_rate - mortal_rate;
            
            // evolution of cover fraction
            m_veget(x,y)[species]+= timestep * growthrate * m_wind_factor;
            
            // limiting conditions
            if(m_veget(x,y)[species] > 1 ){
                m_veget(x,y)[species] = 1;
            }
            if(m_veget(x,y)[species] < 0 || shorefactor == 0){
                m_veget(x,y)[species] = 0;
            }
            
        }
    }
    
    return m_xmin;
}

/* Shift back*/
void vegetation::shiftback(const int plusminus)
{
    m_veget.ShiftOne(plusminus);
}
/*!  Saves the arrays m_u and m_rho.  */
void vegetation::save_arrays(){
    save_2d_vecarray( "veget", m_veget );
}
