/******************************************************************************
  $Id: initsurfgauss.cc,v 1.8 2004/10/01 13:37:28 schatz Exp $
******************************************************************************/

#include "initsurfcone.h"
#include "globals.h"

////////////////////////////////////////
// class CInitSurfCone
//

CInitSurfCone::CInitSurfCone(const dunepar& par, string prefix)
{
  m_volume= par.getrequired<double>(prefix+"cone.volume");
  m_h= par.getrequired<double>(prefix+"cone.h_0");
  m_x0= par.getdefault<double>(prefix+"cone.x_0", 
	      (double)(duneglobals::nx()-1) * duneglobals::dx() / 2.);
  m_y0= par.getdefault<double>(prefix+"cone.y_0", 
		   (double)(duneglobals::ny()-1) * duneglobals::dx() / 2.);

  m_slope = tan(0.9*duneglobals::repose_dyn()*M_PI/180.);
}


/*!  Puts the Cone into the height profile \a array.  */

void CInitSurfCone::init_2d_scal(TFktScal& array)
{
  for (int X=0; X< array.SizeX(); X++)
    for (int Y=0; Y< array.SizeY(); Y++)
      array(X, Y) = 0.0;
    
    double h; 
    if(m_volume!=0) h= pow(3*m_slope*m_slope/M_PI * m_volume, 1./3.);
    else h= m_h;
    cout << "h.cone= " << h << endl;
    const double Ratio = h / m_slope;

    for (int X=0; X< array.SizeX(); X++)
      for (int Y=0; Y< array.SizeY(); Y++){
        const double ratio2= (X*array.Delta() - m_x0)*(X*array.Delta() - m_x0) +
	(Y*array.Delta() - m_y0)*(Y*array.Delta() - m_y0);
	if ( ratio2 <= Ratio*Ratio ){
	    array(X, Y)= h - m_slope*sqrt(ratio2);
	}
	else {
	    array(X, Y)= 0;
	}
      }
}

