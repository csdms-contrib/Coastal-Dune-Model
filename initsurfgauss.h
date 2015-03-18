/******************************************************************************
  $Id: initsurfgauss.h,v 1.4 2005/04/08 11:48:45 schatz Exp $
******************************************************************************/

#ifndef __INITSURFGAUSS_H__
#define __INITSURFGAUSS_H__

#include <vector>

#include "globals.h"
#include "initsurf.h"
#include "func.h"

/*!  Height profile initialisation with a superposition of Gaussian hills.  Any
  number of hills can be given in the parameter file with the parameters
  gauss.h_0, gauss.sigma, gauss.x_0 and gauss.y0.  They are arrays giving the
  height, width and position od the hills.  */

class CInitSurfGauss : public arrayinit
{
public:
  CInitSurfGauss(const dunepar& P, string prefix= "");
  ~CInitSurfGauss() {}

  virtual void init_2d_scal(TFktScal& array);

private:
  double m_h, m_H, m_xslope, m_yslope, m_sigma, m_Xsigma, m_x0, m_y0;
  double m_volume;
  
  /*!  Scale factor for y coordinate to allow creation of gaussian hills with
    elliptical equipotential lines.  Set by parameter gauss.scale_y.  */
  double m_xscale, m_yscale;
  double m_maxslope;
};

#endif
