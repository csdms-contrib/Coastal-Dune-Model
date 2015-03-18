/******************************************************************************
  $Id: initsurfgauss.h,v 1.4 2005/04/08 11:48:45 schatz Exp $
******************************************************************************/

#ifndef __INITSURFWAVES_H__
#define __INITSURFWAVES_H__

#include <vector>

#include "globals.h"
#include "initsurf.h"
#include "func.h"

/*!  Height profile initialisation with a superposition of Gaussian hills.  Any
  number of hills can be given in the parameter file with the parameters
  gauss.h_0, gauss.sigma, gauss.x_0 and gauss.y0.  They are arrays giving the
  height, width and position od the hills.  */

class CInitSurfWaves : public arrayinit
{
public:
  CInitSurfWaves(const dunepar& P, string prefix= "");
  ~CInitSurfWaves() {}

  virtual void init_2d_scal(TFktScal& array);

private:
  double m_A, m_wavelenght, m_exposed, m_asym;
  
};

#endif
