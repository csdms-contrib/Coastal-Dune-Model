/******************************************************************************
 $Id: initsurfgauss.cc,v 1.8 2004/10/01 13:37:28 schatz Exp $
 ******************************************************************************/

#include "initsurfgauss.h"
#include "globals.h"

////////////////////////////////////////
// class CInitSurfGauss
//

CInitSurfGauss::CInitSurfGauss(const dunepar& par, string prefix)
{
    m_h= par.getdefault(prefix+"gauss.h_0", 0.0);
    m_H= par.getdefault(prefix+"gauss.plain_height", 0.0);
    m_Xsigma= par.getdefault(prefix+"gauss.Xsigma", 0);
    m_x0= par.getdefault(prefix+"gauss.x_0",
                         (double)(duneglobals::nx()-1) * duneglobals::dx() / 2.);
    m_y0= par.getdefault(prefix+"gauss.y_0",
                         (double)(duneglobals::ny()-1) * duneglobals::dx() / 2.);
    
    m_xscale = par.getdefault(prefix+"gauss.scale_x", 1.0);
    m_yscale = par.getdefault(prefix+"gauss.scale_y", 1.0);
    m_maxslope = tan(par.getdefault(prefix+"gauss.max.slope", 30.0)*M_PI/180.);
}


/*!  Puts a number of Gaussian hills into the height profile \a array.  */

void CInitSurfGauss::init_2d_scal(TFktScal& array)
{
    for (int iX=0; iX< array.SizeX(); iX++)
        for (int iY=0; iY< array.SizeY(); iY++)
            array(iX, iY) = 0.0;
    
    m_sigma= m_h*exp(-0.5)/m_maxslope;

    for (int iX=0; iX< array.SizeX(); iX++)
        for (int iY=0; iY< array.SizeY(); iY++) {
            const double X= (iX * array.Delta() - m_x0);
            const double dX= (X>0 ? 1. : -1.)*m_Xsigma*array.SizeX()*array.Delta();
            const double Y= (iY * array.Delta() - m_y0);
            double Ratio2= ((X*X/(m_xscale*m_xscale)+Y*Y/(m_yscale*m_yscale)<dX*dX)? 0 : ((X-dX)*(X-dX)/(m_xscale*m_xscale) + (Y-dX)*(Y-dX)/(m_yscale*m_yscale))/(2 *
                                                                                                                                                                  m_sigma*m_sigma));
            
            const double plane= m_H;
            
            array(iX, iY)= m_h * exp( - Ratio2 ) + (plane >= 0 ? plane : 0);
        }
}

