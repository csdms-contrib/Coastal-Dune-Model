/******************************************************************************
  $Id: initsurfmatlab.cc,v 1.4 2004/09/30 10:05:08 schatz Exp $
******************************************************************************/

#include <sstream>

#include "initsurfmatlab.h"

/////////////////////////////////////////////
// Surface from MATLAB file
//

CInitSurfMatlab::CInitSurfMatlab(const dunepar& P, string prefix)
{
  m_filename= P.getrequired<string>(prefix+"matlab.file");
}


void CInitSurfMatlab::Create(PTG_Func2dScalar& f)
{
  ifstream is(m_filename.c_str());
  if (!is) {
    cout << "Opening \"" + m_filename + "\" failed" << endl;
    exit(1);
  }

  // ---- read dimensions, delta ----
  double dLX,dLY,dDelta;
  is >> dLX >> dLY >> dDelta;

  const int iNX = (int)(dLX+0.5);
  const int iNY = (int)(dLY+0.5);

  if (iNX > f.SizeX() || iNY > f.SizeY()) {
    stringstream s;
    s << "The Dimension(s) of the MATLAB file are too large:\n"
      << "Simulation:  (" << f.SizeX() << ", " << f.SizeY() << ")\n"
      << "MATLAB file: (" << iNX << ", " << iNY << ")" << ends;
    cout << s.str() << endl;
    exit(1);
  }

  if (fabs(dDelta-f.Delta()) > 1e-6) {
    stringstream s;
    s << "Discretisation of MATLAB file and simulation is different:\n"
      << "Simulation:  " << f.Delta() << "\n"
      << "MATLAB file: " << dDelta << ends;
    cout << s.str() << endl;
    exit(1);
  }

  // ---- read field ----
  
  const int iN0X = (f.SizeX() - iNX) / 2;
  const int iN0Y = (f.SizeY() - iNY) / 2;

  for (int x=0; x < f.SizeX(); x++) {
    for (int y=0; y < f.SizeY(); y++) {
      if (x >= iN0X && x < iN0X+iNX && y >= iN0Y && y < iN0Y+iNY) {
	is >> f[x][y];
      } else {
	f[x][y] = 0.0;
      }
    }
  }
  if (!is) {
    cout << "Error reading data from: \"" + m_filename + "\"" << endl;
    exit(1);
  }
}






