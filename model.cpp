#include <iostream>
#include <limits>  // For machine epsilon
#include <numeric> // For std::accumulate
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath> // For math functions like pow
#include <algorithm> // For std::max

using namespace std;
// Define a model function with an arbitrary number of parameters
void model(double Vp, std::vector<double> array1, double& RR, double& Ie, double& In, double& Ip, double& Io) {
    // Array input
    double Te = array1[0];
    double Tn = array1[1];
    double Tp = array1[2];
    double np = array1[3];
    double as = array1[4];
    double Vs = array1[5];
    double rss = array1[6];
    double Spap = array1[7];
    double B = array1[8];
    int hh = (int)array1[9];
    double mr = array1[10];
    np = pow(10, np);
    rss = pow(10,rss);
    // Parameter limits
    Te = std::clamp(Te, 0.0, 20.0);
    Tn = std::clamp(Tn, 0.0, 20.0);
    Tp = std::clamp(Tp, 0.0, 20.0);
    np = std::clamp(np, 0.0, 10e20);
    as = std::clamp(as, 0.15, 1e8);
    Vs = std::clamp(Vs, -10.0, 30.0);
    rss = std::clamp(rss, 0.0, 1e-2);
    //std::cout << "----------------------------------------------------------------\n";
    //std::cout << "Te= " << Te << " Tn= " << Tn << " Tp= " << Tp << " np= "  << np   << "\n";
    //std::cout << "np= " << np << " as= " << as << " Vs =" << Vs << " Spap=" << Spap << "\n";
    // Constants
    const double mi = 1.67262192e-27 * mr;
    const double q = 1.60217663e-19;
    const double me = 9.1093837e-31;
    const double pi = 3.14159265359;
    // Probe size
    const double rp = 0.0005;
    const double h = 0.0100;

    // Sheath radius calculation
    double g = Te / Tn;
    double ubp = sqrt(Te / mi *1.8* q * (1 + as) / (1 + g * as));
    double ubn = sqrt(Tp / mi * q);
    double ve = sqrt(8 * Te / pi / me * q);
    //std::cout << "g= " << g << " ubp= " << ubp << " ubn= " << ubn << " ve= "  << ve   << "\n";

//sheath node array
const double eps0=8.854187817e-12;
double rs=0;
double dr=1*1e-7;
int R=floor(1e-3/dr)+100;
//cout << "Rint=" << R << "\n";
double rrp = floor(rp/dr);
double nnn=0;
double u_sheath=0;
double charge_sign=0;

try
{
for(int ggg = 0;ggg<20;ggg++)
{
//cout << "ggg=" << ggg << "\n";
 //computation node

 double u1=0;
 double E1=0;
 double V1=0;
 double u2=0;
 double E2=0;
 double V2=0;
 V1=Vs;
 //charge particle
 if(Vp<Vs)
 {
  u_sheath=ubp;
  u1=ubp;
  charge_sign=1.0;
  nnn=np;
 }
 else
 {
  u_sheath=ubn;
  u1=ubn;
  charge_sign=-1.0;
  nnn=np*as/(1+as);
 }
 int dR=0;
 int rr=0;
 //field calculation (leap-frog)
 E1=E1+q*nnn*u_sheath/eps0/u1*dr*0.5;
 for (int r=1;r<=R*2;r++)
 {
  u2=u1+charge_sign*q*E1/mi/u1*dr;
  if(u2<0){u2=0;}
  E2=E1+((rrp+R)/(rrp+R-r))*charge_sign*nnn*q*u_sheath/eps0/u1*dr-E1/(rrp+R-r);
  V2=V1-E1*dr;
  //cout << " ggg=" << ggg << " r=" << r << " R=" << R << " V=" << V[r+1] << " E=" << E[r+1] << " u=" << u[r+1] << "\n";
  if(Vp+0.1>V2&&V2>Vp-0.1)
  {
    //cout << " ggg=" << ggg << "get r=" << r << " R=" << R << "Vp" << Vp << " " << V[r+1] <<"\n";
    rs=r*dr;
    rr=r;
    break;
  }
  V1=V2;E1=E2;u1=u2;
 }//i

 if(rr<R+2&&rr>R-3)
 {
    break;
 }
 dR = R-rr;
 R = R-floor(dR/1.2);
}
//cout << "try_out" << " rs=" << rs << "R=" << R << "\n";
}

catch (double e)
{
    cout << "An exception occurred. Exception Nr. " << e << '\n';
}
    // Current calculation
    double Se = 0;
    double Sg = 0;
    double b = rp * q * B / sqrt(me * Te * q);
    //std::cout << "R= " << R << " rp= " << rp << " rs= " << rs << "\n";
    //std::cout << "me= " << me << " q= " << q << " B=" << B << " Spap=" << Spap << "\n";
    //std::cout << "b=" << b << "\n";
    Sg = 2 * pi * rp * h + pi * rp * rp;
    Se = Sg * exp(-b * b / 2) + Spap * (1 - exp(-b * b / 2));
    double Sss = 2 * pi * h * (rp + rs) + pi * (rp + rs) * (rp + rs);
    //std::cout << " Sg= " << Sg << " Se= " << Se << " Sss= "  << Sss   << "\n";

    // Current equation
    if (Vp < Vs) {

        //cout << "Low volt\n";
        //cout << "Ie=" <<exp((1 / Te) * (Vp - Vs)) << "\n";
        //cout << "Ie=" << (0.25 * q * ve * np / (1 + as)) << "\n";
        //cout << "Ie=" << Se << "\n";
        Ie = (0.25 * q * ve * np / (1 + as)) * Se * exp((1 / Te) * (Vp - Vs));
        //cout << "Ie=" << Ie << "\n";

        //cout << "Ip=" << Sss << "\n";
        //cout << "Ip=" << (-q * ubp * np) << "\n";
        Ip = (-q * ubp * np) * Sss;
        //cout << "Ip=" << Ip << "\n";

        //cout << "In=" << exp((1 / Tn) * (Vp - Vs)) << "\n";
        //cout << "In=" << (q * ubn * np * as / (1 + as)) << "\n";
        In = (q * ubn * np * as / (1 + as)) * Sg * exp((1 / Tn) * (Vp - Vs));
        //cout << "In=" << In << "\n";

        //cout << "Low Ie=" << Ie << " In=" << In << " Ip=" << Ip << " Io=" << Io << "\n";
    }
    else
    {
        //cout << "High volt\n";
        //cout << "High Ie=" << Ie << " In=" << In << " Ip=" << Ip << " Io=" << Io << "\n";
   //-------------------------------------------------------
        if(hh==1)
        {
          rss=rs;
          //std::cout << "mode 1\n";
          //std::cout << "rss=1 "<< rss << "\n";
        }
        else if(hh==2)
        {
          rss=1;
          //std::cout << "mode 2\n";
          //std::cout << "rss=2 "<< rss << "\n";
        }
        else if(hh==3)
        {
          //std::cout << "mode 3\n";
          //std::cout << "rss=3 "<< rss << "\n";
        }
        else if(hh==4)
        {
          double qq = (1 / Te) * abs(Vp - Vs);
          double a1 = exp(qq);
          double a2 = erfc(sqrt(qq));
          double a3 = 2 * sqrt(qq / 3.14159);
          //std::cout << "mode 4\n";
          //cout << " a1= " << a1 << " a2= " << a2  << " a3= " << a3 << "\n";
          Ie = (0.25 * q * ve * np / (1 + as)) * Se * (a3 + a1 * a2);
          Ip = (-q * ubp * np) * Sg * exp((1 / Tp) * (Vs - Vp));
          In = (q * ubn * np * as / (1 + as)) * Sss;
        }
        else
        {
          std::cout << "Invalid : h" << "\n";
        }
        if(hh<=3)
        {
         double vv=Vp-Vs;
         double uc=pow(abs(vv/Te)*pow(rp,2)/abs(pow(rss+rp,2)-pow(rp,2)),0.5);
         double a1=(rss+rp)/rp*erf(uc);
         double a2=exp(abs(vv/Te));
         double a3=erfc((rss+rp)/rp*uc);
         if(abs(a2)<1e-300||abs(a3)<1e-300){a2=0;a3=0;}
         double a4=a1+a2*a3;
         Ie = (0.25 * q * ve * np / (1 + as)) * Se * a4;
         Ip = (-q * ubp * np) * Sg *exp((1 / Tp) * (Vs - Vp));
         In = (q * ubn * np * as / (1 + as)) * Sss;
        }
        //cout << "High Ie=" << Ie << " In=" << In << " Ip=" << Ip << " Io=" << Io << "\n";
    }
    RR=rs;
    Io = In + Ip + Ie;
    //cout << "fin Ie=" << Ie << "In=" << In << "Ip=" << Ip << "Io=" << Io << "\n";
}
