/******************************************************************************
 $Id: initsurfgauss.cc,v 1.8 2004/10/01 13:37:28 schatz Exp $
 ******************************************************************************/

#include "initsurfwaves.h"
#include "globals.h"

////////////////////////////////////////
// class CInitSurfGauss
//

CInitSurfWaves::CInitSurfWaves(const dunepar& par, string prefix)
{
    m_A= par.getrequired<double>(prefix+"wave.A", 0.0);
    m_wavelenght= par.getdefault<double>(prefix+"wave.lenght", 1.0);
    m_asym= par.getdefault<double>(prefix+"wave.asym", 0.5);
    m_exposed= par.getdefault<double>(prefix+"wave.exposed", 0.2);
//    m_dY= par.getrequired<double>(prefix+"wave.Ydist", 0.1);
}


/*!  Puts a number of Gaussian hills into the height profile \a array.  */

void CInitSurfWaves::init_2d_scal(TFktScal& array)
{
    for (int iX=0; iX< array.SizeX(); iX++)
        for (int iY=0; iY< array.SizeY(); iY++)
            array(iX, iY) = 0.0;
        
    for (int iX=0; iX< array.SizeX(); iX++)
        for (int iY=0; iY< array.SizeY(); iY++) {
            const double X= iX * array.Delta();
            const double Y= iY * array.Delta();
            double Ly = m_wavelenght; // * array.Delta();
            double Lx = Ly/m_asym;
            
            double phase1 = 2*M_PI * (X/Lx + Y/Ly);
            double phase2 = 2*M_PI * (X/Lx - Y/Ly);
            
            // add noise
//            double r = 1+0.1*exp(-((X-Y0)*(X-Y0)+(Y-Y0)*(Y)))
            double Aphase1 = 2*M_PI * (X/Lx + Y/Ly)/3.;
            double Aphase2 = 2*M_PI * (X/Lx - Y/Ly)/3.;
            double Ar = 0.25*((cos(Aphase1)+cos(Aphase2)));
            
            
            double value = 0.5*(2*(1-m_exposed)-(cos(phase1)+cos(phase2)))/(2-m_exposed);
            array(iX, iY)= m_A * (1 + 0.5*Ar) * (value > 0 ? value : 0);
//            array(iX, iY)= m_A * 0.5*(1-cos(phase));
        }
}

