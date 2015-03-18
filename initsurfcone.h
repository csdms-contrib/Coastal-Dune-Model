/******************************************************************************
  $Id: initsurfgauss.h,v 1.3 2004/08/18 10:58:09 schatz Exp $
******************************************************************************/

#ifndef __INITSURFCONE_H__
#define __INITSURFCONE_H__

#include "globals.h"
#include "initsurf.h"
#include "func.h"


/*!  Height profile initialisation with a superposition of Gaussian hills.  Any
  number of hills can be given in the parameter file with the parameters
  cone.h_0, cone.x_0 and cone.y0.  They is an array giving the
  height, width and position of the cone.  */

class CInitSurfCone : public arrayinit
{
public:
  CInitSurfCone(const dunepar& P, string prefix= "");
  ~CInitSurfCone() {}

  virtual void init_2d_scal(TFktScal& array);

private:
  double m_volume, m_h, m_x0, m_y0;
  
  double m_slope;
};

#endif
