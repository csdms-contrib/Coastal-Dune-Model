#ifndef __INITSURFMATLAB_H__
#define __INITSURFMATLAB_H__

#include "globals.h"
#include "initsurf.h"


// ---- Surface from MATLAB file ----

class CInitSurfMatlab : public arrayinit
{

public:
  CInitSurfMatlab(const dunepar& P, string prefix= "");

  virtual void init_2d_scal(TFktScal& array)
	{ Create(array); }

  virtual void Create(PTG_Func2dScalar& f);

private:
  string m_filename;
};

#endif
