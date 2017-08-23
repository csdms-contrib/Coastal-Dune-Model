/******************************************************************************
 $Id: vegetation.h,v 1.5 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#ifndef VEGETATION_H
#define VEGETATION_H

#include "globals.h"
#include "func.h"
//#include "PTG_FileNames.h"
#include "vec.h"

class arrayinit;

class vegetation : public dunedata
{

public:
    // construction
    vegetation(const dunepar& par);
    virtual ~vegetation() {}
        
    //  Functions reimplemented from coast:
    virtual void init(const dunepar& par);
    virtual int evol(TFktScal& rho_veget, const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dhdt, TFktScal& overwash);
    void getcover(TFktScal& rho_veget);

    virtual void shiftback(const int plusminus);

    virtual void save_arrays();

private:
    
    virtual int evolspec(const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dhdt, TFktScal& overwash, int species);

    /* auxiliary function for vegetation */
    TFktVec m_veget;    // vegetation cover fraction for two species

    /*! Vegetation parameters (see default.par)*/
    int m_xmin, m_xmin0; // vegetation minimum distance to shoreline
    int m_veg_type;
    double m_Lveg; // vegetation limit = m_xmin * dx
    double m_zmin; // veget minimum height above sealevel
    double m_veget_init0, m_veget_init1;
    double m_Tveg, m_rho_max, m_rho_min, m_Hveg, m_sens, m_veg_season, m_wind_factor, m_Vlateral_factor, m_angle_ref;
    bool m_spec1, m_spec2, m_lateral;
    bool m_survive;

};


#endif //  VEGETATION_H

