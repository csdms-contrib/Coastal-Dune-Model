/******************************************************************************
  $Id: initsurfparabol.h,v 1.4 2004/09/30 10:05:09 schatz Exp $
******************************************************************************/


#ifndef __INITSURFPARABOL_H__
#define __INITSURFPARABOL_H__

#include "globals.h"
#include "initsurf.h"

// ---- Paraboloid / Hyperboloid Barchan ----

class CInitSurfParabol : public arrayinit
{

public:
  CInitSurfParabol(const dunepar& P, string prefix= "");

  void CreateSlipFace(int iNX, int iNY, double dDelta);
  PTG_Func<double>& GetSlipFace() { return m_fSlipFace; }

  virtual void init_2d_scal(TFktScal& array)
	{ Create(array); }

  virtual void Create(PTG_Func2dScalar& f);

private:
  PTG_Func<double> m_fSlipFace;
  bool m_bSlipFaceCreated;

  double dTheta;
  double dH;

  double dL;
  double dW;

  double dB0;
  double dB2;
  
  double m_x0, m_y0;
  bool m_cos, m_sym;
};

#endif
