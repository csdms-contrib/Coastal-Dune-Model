#ifndef __ANALYZE_NEW_H__
#define __ANALYZE_NEW_H__

#include "globals.h"
#include "func.h"


class analyze
{
public:
  // creation
  analyze(const dunepar& P);

  void Calc(int t, double timestep,
  	double qin, double qout, const TFktScal& h, const TFktScal& m_rhoveg);

  double Center() { return m_dCenter; }

private:
  // internal helper functions
  double GetMaxPos(double f0, double f1, double f2);

private:
  double m_dCenter;
  double m_zmin;

  ofstream m_os;
};

#endif
