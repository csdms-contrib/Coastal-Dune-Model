#ifndef __ANALYZE_H__
#define __ANALYZE_H__

#include "globals.h"
#include "func.h"


class CAnalyze
{
public:
  // creation
  CAnalyze(const dunepar& P,
	   const TFktScal& h, const TFktScal& dh, const TFktScal& ex, const TFktVec& flux,
	   const TFktVec& tau);

  void Calc(int t, double dDT, double dInFlux, double dOutFlux,
	    int iShiftX, double dVDune, 
	    double MeanRho, int iIterRho, int iIterAval);

  double Center() { return m_dCenter; }

private:
  // internal helper functions
  double GetMaxPos(double f0, double f1, double f2);

private:
  const dunepar& m_P;

  const TFktScal& m_h;
  const TFktScal& m_dh;
  const TFktScal& m_ex;
  const TFktVec& m_flux;
  const TFktVec& m_tau;

  double m_dTime;
  double m_dCenter;

  bool m_bMsg;
  ofstream m_os;

  double m_dRhoSand;
};

#endif
