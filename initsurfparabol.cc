/******************************************************************************
  $Id: initsurfparabol.cc,v 1.4 2004/09/30 10:05:08 schatz Exp $
******************************************************************************/

#include <sstream>

#include "initsurfparabol.h"

/*
Create a barchan dune as a parabola.
*/

/////////////////////////////////////////////
// Paraboloid / Hyperboloid Parabol
//

CInitSurfParabol::CInitSurfParabol(const dunepar& P, string prefix) :
  m_bSlipFaceCreated(false)
{
  dTheta  = P.getdefault(prefix+"parabola.theta", 31.0);
  dH = P.getdefault(prefix+"parabola.h", 5.0);
  dL= P.getdefault(prefix+"parabola.l", 5.4*dH + 10.7 );
  dW= P.getdefault(prefix+"parabola.w", 5.8*dH + 1.9 );
  dB0= P.getdefault(prefix+"parabola.b0", -2.4*dH + 7.9 );
  dB2= P.getdefault(prefix+"parabola.b2", 1.46 * dL / (dW*dW) );
  m_x0 = P.getdefault(prefix+"parabola.x_0", duneglobals::dx()*(duneglobals::nx()-1)/2.);
  m_y0 = P.getdefault(prefix+"parabola.y_0", duneglobals::dx()*(duneglobals::ny()-1)/2.);
  m_cos = P.getdefault(prefix+"parabola.cos", false);
  m_sym= P.getdefault(prefix+"parabola.sym", false);
}


void CInitSurfParabol::CreateSlipFace(int iNX, int iNY, double dDelta)
{
  if(m_bSlipFaceCreated) {
    cerr << "CInitSurfParabol::CreateSlipFace: slip face was created before?!\n";
    return;
  }
  // ---- calc slip face ----

  m_fSlipFace.Create(iNY, dDelta);
  for (int y=0; y < iNY; y++) {
    m_fSlipFace[y]= m_x0 + dB0 + dB2 * (y * dDelta - m_y0)*(y * dDelta - m_y0);
  }

  m_bSlipFaceCreated = true;
}



void CInitSurfParabol::Create(PTG_Func2dScalar& f)
{
  const int iNX = f.SizeX();
  const int iNY = f.SizeY();
  const double dDelta = f.Delta();

  if (!m_bSlipFaceCreated) {
    CreateSlipFace(iNX, iNY, dDelta);
  }

  if (m_cos) {
    cout << "Creating cosine surface ..." << endl;
  } else {
    cout << "Creating parabolic surface ..." << endl;
  }

  // ---- calc surface ----

  for (int x=0; x < iNX; x++) {
    for (int y=0; y<iNY; y++) {
      const double dX2 = (x*dDelta - m_x0 )*(x*dDelta - m_x0 ) / (dL*dL);
      const double dY2 = (y*dDelta - m_y0 )*(y*dDelta - m_y0 ) / (dW*dW);
      if (m_cos) {
	// cos: z = h*( cos(sqrt(2)*r) ) )
        //      r = sqrt( (x/l)^2 - (y/w)^2 )
	f(x,y)= dH * cos( sqrt(2*(dX2+dY2)) );
      }
      else {
	// parabol: z = h*(1 - (x/l)^2 - (y/w)^2)
	f(x,y) = dH * (1 - dX2 - dY2);
      }
      if( f(x, y)<= 0.0 )
	f(x, y)= 0.0;
    }
  }


  // ---- equilibrate -> create slip face by removing mass ----
  
  const double dTanTheta = tan(dTheta * M_PI/180.0);
  const double dTanThetaLimit = tan((dTheta+1) * M_PI/180.0);

  PTG_Func2dVec grad(iNX, iNY, dDelta);
  PTG_Func2dVec grad_mid(iNX, iNY, dDelta);
  PTG_Func2dScalar grad_abs(iNX, iNY, dDelta);
  grad_abs.SetAll(0);

  // iteration loop for the mass removement
  for (int i=0; i<1000; i++) {
    grad_mid.GradMid( f );

    for (int x=1; x<iNX-1; x++) {
      for (int y=1; y<iNY-1; y++) {
	grad_mid(x,y) = -grad_mid(x,y);
      }
    }

    grad.GradUpWind( f, grad_mid );

    for (int x=1; x<iNX-1; x++) {
      for (int y=1; y<iNY-1; y++) {
	grad_abs(x,y) = sqrt( grad(x,y)[0]*grad(x,y)[0] + 
				grad(x,y)[1]*grad(x,y)[1]);
      }
    }

    // diagnostic / debug
    if (grad_abs.GetMax() > dTanThetaLimit) {
      stringstream s;
      s << "warning: over relaxation: max(grad)="
	<< grad_abs.GetMax() << ends;
      //      throw CError(s.str());
      cout << s.str() << endl;
      exit(1);
    }

    // change function
    double dChange = 0.0;
    for (int x=1; x<iNX-1; x++) {
      for (int y=1; y<iNY-1; y++) 
	{
	  // A parabola defines the crest, so mass is only
	  // removed in the area inside of the parabola.
	  if ( x*dDelta > m_fSlipFace.GetInterpolatedValue(y*dDelta) )
	    {
	      // Mass is only removed until the angle of repose 
	      // is reached.
	      if (grad_abs(x,y) < dTanTheta) {
		double dz = 0.4 * (dTanTheta-grad_abs(x,y)) * dDelta; 

		// Remove mass only until the zero level is reached.
		if (dz > f(x,y)) {
		  if( dChange < f(x, y) )
		    dChange = f(x,y);
		  f(x,y) = 0.0;
		} else {
		  if( dChange < dz )
		    dChange = dz;
		  f(x,y) = f(x,y) - dz;
		}
	      } // gradient
	    } // parabola
	} // loop
    }
    cout << "i=" << i << ", dz=" << dChange 
	 << ", max(grad(x,y))=" << grad_abs.GetMax() << endl;

    if (dChange < 0.00001) {
      // solution is converged
      break;
    }
  }


  if ( m_sym ) {
    for (int y=0; y < f.SizeY()/2; y++) {
      for (int x=0; x < f.SizeX(); x++) {
	f[x][y] = f[x][f.SizeY()-1-y];
      }
    }
  }

}







