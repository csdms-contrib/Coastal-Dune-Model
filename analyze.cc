#include <string>

#include "analyze.h"

using std::string;

//////////////////////////////////////
// helper classes

struct CMinMax
{
  int m_iMin;
  int m_iMax;
  double m_dMin;
  double m_dMax;
  double m_dh;

  CMinMax(int iMax, double dh) :
    m_iMin(iMax), m_iMax(0),
    m_dMin(iMax), m_dMax(0),
    m_dh(dh)
  {}

  
  double Interpol(double h0, double h1, double href) {
    if (h1-h0 > 1e-10) {
      double retval;
      retval= (href-h0)/(h1-h0);
      if( retval < 0 )
        retval= 0;
      return (retval> 1.0? 1.0 : retval);
    }
    else
      return 0.5;
  }
  

  void operator()(int i, double dl, double d, double dr)
  {
    if (i<=m_iMin) {
      m_iMin = i;
      if( m_dMin > i-1+Interpol(dl,d,m_dh) )
        m_dMin = i-1+Interpol(dl,d,m_dh);
    }

    if (i>=m_iMax) {
      m_iMax = i;
      if( m_iMax < i+Interpol(d,dr,m_dh) )
        m_dMax = i+Interpol(d,dr,m_dh);
    }
  }

  void Check(const char* pszName)
  {
    if (m_dMin > m_dMax) {
      double d = m_dMin;
      m_dMin = m_dMax;
      m_dMax = d;
      cout << pszName << ": Min/max detection failed!" << endl;
    }
  }

  int Center() {
    return (m_iMin+m_iMax)/2;
  }

};





//////////////////////////////////////
// analyze
//

CAnalyze::CAnalyze(const dunepar& P,
		   const TFktScal& h, 
		   const TFktScal& dh,
		   const TFktScal& ex,
		   const TFktVec& flux,
		   const TFktVec& tau)
:
  m_P(P),
  m_h(h),
  m_dh(dh),
  m_ex(ex),
  m_flux(flux),
  m_tau(tau),
  m_dTime(0.),
  m_bMsg(true)
{
  m_dRhoSand = m_P.getdefault("rho_sand", 1650);

  std::string strFile = m_P.getdefault<string>("save.dir", "./");
  if (strFile[strFile.length()-1] != '/') {
    strFile += '/';
  }
  strFile += "time.dat";
  m_os.open(strFile.c_str());

  if (!m_os) {
    cout << "Open file \"" << strFile << "\" failed!" << endl;
    exit(1);
  }

  m_os << "# 1: iteration n\n"
       << "# 2: time t\n"
       << "# 3: height h\n"
       << "# 4: width w\n"
       << "# 5: length l\n"
       << "# 6: horn length l_h\n"
       << "# 7: slip-face length l_s\n"
       << "# 8: windward side length l_0\n"
       << "# 9: volume V\n"
       << "#10: crest position x_c\n"
       << "#11: flux at the crest q_x\n"
       << "#12: center of mass x_m\n"
       << "#13: number of shift backs\n"
       << "#14: v_dune\n"
       << "#15: dune in flux\n"
       << "#16: dune out flux\n"
       << "#17: erosion-deposition integral\n"
       << "#18: change in the surface\n" 
       << "#19: N_iter_u_s\n"
       << "#20: N_iter_rho_s\n"
       << "#21: N_iter_aval\n"
       << "#22: time step for dh/dt\n"
       << "#23: pos. of max. h\n"
       << "#24: pos. of max. tau_x\n"
       << "#25: pos. of max. q\n";
  // two-dimensional model
  if (m_P.getdefault("NY", 100) == 3) {
    m_os << "#26: mean dist between crests\n"
	 << "#27: number of crests (max. of h)" << endl;
  } else {
    m_os   << endl;
  }
}


void CAnalyze::Calc(int t, double dDT, double dInFlux, double dOutFlux,
		    int iShiftX, double dVDune, 
		    double MeanRho, int iIterRho, int iIterAval) 
{
  const double dx = m_h.Delta();

  m_dTime += dDT;

  // calc dune properties

  const double dHMax = m_h.GetMax();
  const double dHRef = dHMax/100.;
  
  double minPos = 0;
  for (int x=m_h.SizeX()-1; x > 0; x--) 
  	if(m_h(x,m_h.SizeY()/2) > 0.1){
		minPos = x;
		break;
	}

  CMinMax mmX(m_h.SizeX(), dHRef);
  CMinMax mmY(m_h.SizeY(), dHRef);

  // scan surface at level dHRef -> length and width

  for (int y=1; y<m_h.SizeY()-1; y++) {
    for (int x=1; x<m_h.SizeX()-1; x++) {
      if (m_h(x,y) > dHRef) {
	mmX(x, m_h(x-1,y), m_h(x,y), m_h(x+1,y));
	mmY(y, m_h(x,y-1), m_h(x,y), m_h(x,y+1));
      }
    }
  }
  mmX.Check("Length L");
  mmY.Check("Width W");


  // scan symmetry line -> l_h, l_0

  CMinMax mmXSym(m_h.SizeX(), dHRef);

  const int iYSym = mmY.Center();
  for (int x=1; x<m_h.SizeX()-1; x++) {
    if (m_h(x,iYSym) > dHRef) {
      mmXSym(x, m_h(x-1,iYSym), m_h(x,iYSym), m_h(x+1,iYSym));
    }
  }
  mmXSym.Check("Symmetry line L0+Ls");


  // calc position of max height, flux, shear

  double dMaxH = 0.;
  double dMaxQ = 0.;
  double dMaxTau = 0.;
  int iXH(0), iXQ(0), iXTau(0);
  for (int x=mmXSym.m_iMin+1; x<mmXSym.m_iMax-2; x++) 
  {
    if (dMaxH < m_h(x,iYSym)) {
      dMaxH = m_h(x,iYSym);
      iXH = x;
    }

    if (dMaxTau < m_tau(x,iYSym)[0]) {
      dMaxTau = m_tau(x,iYSym)[0];
      iXTau = x;
    }

    // upwind
    double dQ = m_flux(x,iYSym)[0];
    if (dMaxQ < dQ) {
      dMaxQ = dQ;
      iXQ = x;
    }
  }
  double dXH = GetMaxPos(m_h(iXH-1,iYSym), m_h(iXH,iYSym), m_h(iXH+1,iYSym));
  double dXTau = GetMaxPos(m_tau(iXTau-1,iYSym)[0], m_tau(iXTau,iYSym)[0], m_tau(iXTau+1,iYSym)[0]);
  double dXQ = GetMaxPos(m_flux(iXQ-1,iYSym)[0],
			 m_flux(iXQ,iYSym)[0],
			 m_flux(iXQ+1,iYSym)[0]);

  dXH = ((double)iXH + dXH) * dx;
  dXQ = ((double)iXQ + dXQ) * dx + 0.5*dx;        // q has been calculated on the faces
  dXTau = ((double)iXTau + dXTau) * dx - 0.5*dx;  // tau_x has been calculated on the faces

  // scan for the crest position

  int iXBrink = -1;
  double m_sepslope = tan(M_PI*m_P.getdefault("sep.angle", 30.0)/180.0);
  for (int x=mmXSym.m_iMin; x<mmXSym.m_iMax-1; x++) {
    if (m_h(x,iYSym)-m_h(x-1,iYSym) > -m_sepslope*dx && m_h(x+1,iYSym)-m_h(x,iYSym) < -m_sepslope*dx) {
      iXBrink = x;
      break;
    }
  }

  if (iXBrink == -1) {
    // use maximum instead / no slip-face
    double dMax = 0.;
    for (int x=mmXSym.m_iMin; x<mmXSym.m_iMax-1; x++) {
      if (dMax < m_h(x,iYSym)) {
	dMax = m_h(x,iYSym);
	iXBrink = x;
      }
    }
  }

  if (iXBrink == -1) {
    iXBrink = (mmXSym.m_iMin+mmXSym.m_iMax)/2;
    if (m_bMsg) {
      cout << "--> Brink/crest detection failed!" << endl;
    }
  }

  // ---- calc lengths ----

  m_dCenter = (mmX.m_dMax+mmX.m_dMin)*dx/2;

  const double dL = (mmX.m_dMax-mmX.m_dMin)*dx;
  const double dW = (mmY.m_dMax-mmY.m_dMin)*dx;
  const double dL0 = (iXBrink-mmXSym.m_dMin)*dx;
  const double dLs = (mmXSym.m_dMax-iXBrink)*dx;
  const double dLh = (mmX.m_dMax-mmXSym.m_dMax)*dx;


  // ---- calc fluxes ----

  const double dFluxBrink = m_flux(iXBrink,iYSym)[0];
  
  double MeanFlux = 0;
  int N = 0;
  for (int y=0; y<m_h.SizeY(); y++) {
    for (int x=0; x<m_h.SizeX(); x++) {
    	if(m_h(x,y)>1e-2){
    		MeanFlux += sqrt(m_flux(x,y)[0]*m_flux(x,y)[0] + m_flux(x,y)[1]*m_flux(x,y)[1]);
		N++;
	}
    }
  }
  MeanFlux /= N;
  
  // ---- volume / mass ----

  const double dVol = m_h.Integrate();
  
  const double dVoldt = m_dh.Integrate()/dDT;
  
  // ---- erosion-deposition integral ----
  
  double dVol_ex = m_ex.Integrate(); dVol_ex/= dx*m_h.SizeY();

  // ---- for 1dim model calc crest mean crest distance ----
  double dMeanDist = 0.0;
  int LastPos=0;
  int count=0;
  if (m_h.SizeY() == 3) {
    for ( int i=3; i<m_h.SizeX(); i++) {
      if (m_h(i,1)-m_h(i-1,1) <= 0.0 && m_h(i-1,1)-m_h(i-2,1) > 0.0) {
	count++;
	if (count > 1) {
	  dMeanDist += ((double) (i-LastPos)) * dx; 
	  LastPos = i;
	}
      }
    }
    if (count > 1) dMeanDist /= (double) (count-1);
  }
  // ---- write data ----

  m_os << t << " "                    // 1: iterations
       << m_dTime << " "              // 2: time in s
       << dHMax << " "                // 3: height
       << dW << " "                   // 4: width w
       << dL << " "                   // 5: length l
       << dLh << " "                  // 6: horn length l_h
       << dLs << " "                  // 7: slip-face length l_s
       << dL0 << " "                  // 8: windward side length l_0
       << dVol << " "                 // 9: volume / mass of sand
       << iXBrink*dx << " "           //10: brink position x_c
       << dFluxBrink << " "           //11: flux at the brink q_x
       << m_h.CenterOfMassX() << " "  //12: center of mass
       << iShiftX << " "              //13: number of shift backs
       << dVDune << " "               //14: v_dune
       << dInFlux << " "              //15: dune in flux
       << dOutFlux << " "             //16: dune out flux
       << dVol_ex << " "	      //17: erosion-deposition integral
       << dVoldt << " "		      //18: change in the surface
       << 0 << " "               //19: N_iter_u_s
       << iIterRho << " "             //20: N_iter_rho_s
       << iIterAval << " "            //21: N_iter_aval
       << dDT << " "                  //22: time step for dh/dt
       //<< dXH << " "                  //23: pos. of max. h
       << minPos << " "                  //23: pos. of max. h
       << MeanRho << " "
       << MeanFlux << " "
       << dXTau << " "                //24: pos. of max. tau_x
       << dXQ << " ";                  //25: pos. of max. q
  // two-dimensional model
  if(m_h.SizeY() == 3) {
    m_os << " " << dMeanDist << " "        //26: mean dist. between crests
	 << count << endl;                 //27: number of crests (max. of h)
  }
  else {
    m_os << endl;
  }
}


double CAnalyze::GetMaxPos(double f0, double f1, double f2) 
{
  // assert that: f0 < f1 and f1 > f2
  if (f0 > f1 || f2 > f1)
    return 0;

  // calc derivatives: d1 > 0, d2 < 0
  double d1 = f1-f0;
  double d2 = f2-f1;

  // calc zero pos.  d1 + epsilon * (d2-d1) = 0
  double epsilon;
  if (d1 != d2) {
    epsilon = d1 / (d1-d2);
  } else {
    epsilon = 0.5;
  }
  return epsilon - 0.5;
}
