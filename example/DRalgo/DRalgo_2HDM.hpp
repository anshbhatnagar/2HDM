// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef POTENTIAL_DRALGO_2HDM_HPP_INCLUDED
#define POTENTIAL_DRALGO_2HDM_HPP_INCLUDED

/**
   The xSMin  DRalgo
   https://arxiv.org/pdf/22xx.xxxxx.pdf

  Sqrt -> sqrt, pow -> pow, 2 -> 2.
  M_PI -> M_PI, log-> log
   /.{\[CurlyPhi]^2 -> Hsq, \[CurlyPhi]^4 -> Hsq^2, s^2->Ssq, s^4->Ssq^2,lambdaH->lamH,lambdaHS->lamHS, lambdaS->lamS,\[Mu]Ssq->muSsq,\[Mu]Hsq->muHsq,g1^2 -> g1sq, g2^2 -> g2sq, g13d^2->g13dsq,g23d^2->g23dsq,g33d^2->g33dsq,lambdaH3d->lamH3d, lambdaHS3d->lamHS3d,lambdaS3d->lamS3d}
 

*/

#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Eigenvalues>
#include <interpolation.h>

#include "pow.hpp"
#include "potential.hpp"
#include "DRalgo_2HDM_RenormParam.hpp"
#include "DRalgo_2HDM_Potential.hpp"
#include "models/SM_parameters.hpp"

namespace EffectivePotential {

class DR_2HDM: public Potential {
  public :

    DR_2HDM(std::vector<double> input) {
    
      bool print = false;

      double mZ = 91.1876;

      DR_2HDM_RenormParam params(input);

      InputParams = input;
      
      /*if(print){
        std::cout << "Input Ms (physical scalar mass) = " << Ms << std::endl;
        std::cout << "Input lamHS = " << lamHS_input << std::endl;
        std::cout << "Input lamS = " << lamS_input << std::endl;
      }
      
      if(print){
        std::cout << "msSq = " << msSq_input << std::endl;
        std::cout << "mphiSq = " << mphiSq_input << std::endl;
        std::cout << "lamH = " << lamH_input << std::endl;
      }*/

      // Check points
      /*if( !check_points(print) ) { 
        std::cout << "Failed!" << std::endl;
        exit(EXIT_FAILURE);
      }*/

      LagrangianParams = params.OneLoopRenormalised();
  
      solveBetas(LagrangianParams, input[7], 1., 5000.);

    }

    size_t get_n_scalars() const override { return 2;}

    void printLagrangianParams(){
      std::vector<std::string> paramNames = {"g1Sq", "g2Sq", "g3Sq", "l1", "l2", "l3", "l4", "l5", "sqrt(ytSq)", "m12RSq", "m1Sq", "m2Sq"};

      for(int i = 0; i < LagrangianParams.size(); i++){
        std::cout << paramNames[i] << ": " << LagrangianParams[i] << std::endl;
      } 
    }

    void printHiggsvev(){
      double l4 = LagrangianParams[6];
      double l5 = LagrangianParams[7];
      double lF = l5 - l4;

      double mHpm = InputParams[3];
      double mA = InputParams[2];

      std::cout << "Higgs vev is: " << sqrt(2.*(pow(mHpm,2.) - pow(mA,2.))/lF) << std::endl;
    }

    std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
      auto phi1 = phi;
      phi1[0] = - phi[0];
      auto phi2 = phi;
      phi2[1] = - phi[1];
      return {phi1,phi2};
    };

    /*bool check_points(bool print_) {

      if ( mphiSq_input > 0){
        std::cout << "Higgs mass (mphiSq) is positive." << std::endl;
        return false;
      }

      if ( msSq_input > 0){
        std::cout << "Scalar mass (msSq) is positive." << std::endl;
        return false;
      }

      if(print_) {
        std::cout << "- pow(mphiSq_input,2)/lamH_input = " << - pow(mphiSq_input,2)/lamH_input << std::endl;
        std::cout << "- pow(msSq_input,2)/lamS_input = " << - pow(msSq_input,2)/lamS_input << std::endl;
      }

      return true;
    }*/

    /*double get_TC_from_expression() const {
      const double cs = 1. / 12. * (2. * lamHS_input + 3. * lamS_input);
      const double ch = 1. / 48. * (9. * g2sq_input + 3. * g1sq_input + 12. * yt_input*yt_input + 24. * lamH_input + 2. * lamHS_input);
      const double TC_sq = -(lamS_input * ch * mphiSq_input - lamH_input * cs * msSq_input + std::sqrt(lamS_input * lamH_input) * std::abs(cs * mphiSq_input - ch * msSq_input)) / (lamS_input * square(ch) - lamH_input * square(cs));
      return std::sqrt(TC_sq);
      }*/

    double V(Eigen::VectorXd phi, double T) const override {
    
      const std::vector<std::complex<double>> par = get_3d_parameters(T);

      DR_2HDM_Potential V2HDM;
      std::complex<double> veffLO = V2HDM.V0(phi, T, par);
      std::complex<double> veffNLO = V2HDM.V1(phi, T, par);
      std::complex<double> veffNNLO = V2HDM.V2(phi, T, par);

      double Veff = veffLO.real() + veffNLO.real() + veffNNLO.real();

      return T * Veff;
      
    }

    void Betas( const std::vector<double>& x, std::vector<double>& dxdt, const double t) override {
    
      double g1Sq = x[0];
      double g2Sq = x[1];
      double g3Sq = x[2];
      double lam1H = x[3];
      double lam2H = x[4];
      double lam3H = x[5];
      double lam4H = x[6];
      double lam5H = x[7];
      double yt2 = x[8];
      double m12RSq = x[9];
      double m1Sq = x[10];
      double m2Sq = x[11];
      dxdt[0] = 1./t * (7*pow(g1Sq,2.))/(8.*pow(pi,2));
      dxdt[1] = 1./t * (-3*pow(g2Sq,2.))/(8.*pow(pi,2));
      dxdt[2] = 1./t * (-7*pow(g3Sq,2.))/(8.*pow(pi,2));
      dxdt[3] = 1./t * (3*pow(g1Sq,2.) + 9*pow(g2Sq,2.) + 6*g1Sq*(g2Sq - 4*lam1H) - 72*g2Sq*lam1H + 8*(24*pow(lam1H,2) + 2*pow(lam3H,2) + 2*lam3H*lam4H + pow(lam4H,2) + pow(lam5H,2)))/(128.*pow(pi,2));
      dxdt[4] = 1./t * (3*pow(g1Sq,2.) + 9*pow(g2Sq,2.) + 6*g1Sq*(g2Sq - 4*lam2H) - 72*g2Sq*lam2H + 96*lam2H*(pow(yt2,2) + 2*lam2H) + 8*(-6*pow(yt2,4) + 2*pow(lam3H,2) + 2*lam3H*lam4H + pow(lam4H,2) + pow(lam5H,2)))/(128.*pow(pi,2));
      dxdt[5] = 1./t * (3*pow(g1Sq,2.) + 9*pow(g2Sq,2.) - 36*g2Sq*lam3H - 6*g1Sq*(g2Sq + 2*lam3H) + 8*(lam3H*(3*pow(yt2,2) + 6*(lam1H + lam2H) + 2*lam3H) + 2*(lam1H + lam2H)*lam4H + pow(lam4H,2) + pow(lam5H,2)))/(64.*pow(pi,2));
      dxdt[6] = 1./t * (3*g1Sq*(g2Sq - lam4H) - 9*g2Sq*lam4H + 6*pow(yt2,2)*lam4H + 4*lam4H*(lam1H + lam2H + 2*lam3H + lam4H) + 8*pow(lam5H,2))/(16.*pow(pi,2));
      dxdt[7] = 1./t * ((-3*g1Sq - 9*g2Sq + 6*pow(yt2,2) + 4*(lam1H + lam2H + 2*lam3H + 3*lam4H))*lam5H)/(16.*pow(pi,2));
      dxdt[8] = 1./t * (yt2*(-17*g1Sq - 27*g2Sq - 96*g3Sq + 54*pow(yt2,2)))/(192.*pow(pi,2));
      dxdt[9] = 1./t * (m12RSq*(-3*g1Sq - 9*g2Sq + 6*pow(yt2,2) + 4*lam3H + 8*lam4H + 12*lam5H))/(32.*pow(pi,2));
      dxdt[10] = 1./t * (-3*m1Sq*(g1Sq + 3*g2Sq - 8*lam1H) + 4*m2Sq*(2*lam3H + lam4H))/(32.*pow(pi,2));
      dxdt[11] = 1./t * (-3*m2Sq*(g1Sq + 3*g2Sq - 4*(pow(yt2,2) + 2*lam2H)) + 4*m1Sq*(2*lam3H + lam4H))/(32.*pow(pi,2));
    
    }

    std::vector< std::complex<double>> get_3d_parameters(double T) const {
      double mu4D = 1;
      double mu3D = 1;
      //double Gamma = scaleFactor * T;
      //double Lb = 2.*log(scaleFactor) + 2. * EulerGamma - 2. * log(4 * M_PI);
      //double Lf = Lb + 4. * log(2.);
      double EulerGamma = 0.5772156649;
      double Glaisher = 1.282427;
      std::complex<double> mu4 = 4 * pi * exp(EulerGamma) * mu4D * T;
      std::complex<double> mu3 = mu3D * T;
      std::complex<double> Lb = 2*EulerGamma - 2*log(4*pi) + log(pow(mu4,2)/pow(T,2));
      std::complex<double> Lf = 2*EulerGamma + 4*log(2) - 2*log(4*pi) + log(pow(mu4,2)/pow(T,2));
      //std::complex<double> g1Sq, g2Sq, g3Sq, lam1H, lam2H, lam3H, lam4H, lam5H, yt2, m12RSq, m1Sq, m2Sq;
      
      std::complex<double> g1Sq(alglib::spline1dcalc(RGEs[0], mu4.real()), 0.);
      std::complex<double> g2Sq(alglib::spline1dcalc(RGEs[1], mu4.real()), 0.);
      std::complex<double> g3Sq(alglib::spline1dcalc(RGEs[2], mu4.real()), 0.);
      std::complex<double> lam1H(alglib::spline1dcalc(RGEs[3], mu4.real()), 0.);
      std::complex<double> lam2H(alglib::spline1dcalc(RGEs[4], mu4.real()), 0.);
      std::complex<double> lam3H(alglib::spline1dcalc(RGEs[5], mu4.real()), 0.);
      std::complex<double> lam4H(alglib::spline1dcalc(RGEs[6], mu4.real()), 0.);
      std::complex<double> lam5H(alglib::spline1dcalc(RGEs[7], mu4.real()), 0.);
      std::complex<double> yt2(alglib::spline1dcalc(RGEs[8], mu4.real()), 0.);
      std::complex<double> m12RSq(alglib::spline1dcalc(RGEs[9], mu4.real()), 0.);
      std::complex<double> m1Sq(alglib::spline1dcalc(RGEs[10], mu4.real()), 0.);
      std::complex<double> m2Sq(alglib::spline1dcalc(RGEs[11], mu4.real()), 0.);
      
      std::complex<double> g1 = sqrt(g1Sq);
      std::complex<double> g2 = sqrt(g2Sq);
      std::complex<double> g3 = sqrt(g3Sq);

      std::complex<double> m1 = m1Sq;
      std::complex<double> m2 = m2Sq;
      std::complex<double> m12R = m12RSq;

      //---------------------------------------------------------------------------------
      std::complex<double> g13d2 = pow(g1,2)*T - (pow(g1,4)*(Lb + 20*Lf)*T)/(48.*pow(pi,2));
      std::complex<double> g23d2 = pow(g2,2)*T + (pow(g2,4)*(2 + 21*Lb - 12*Lf)*T)/(48.*pow(pi,2));
      std::complex<double> g33d2 = pow(g3,2)*T + (pow(g3,4)*(1 + 11*Lb - 4*Lf)*T)/(16.*pow(pi,2));
      std::complex<double> lam1H3d = (T*((pow(g1,4) + 2*pow(g1,2)*pow(g2,2) + 3*pow(g2,4))*(2 - 3*Lb) + 24*(pow(g1,2) + 3*pow(g2,2))*Lb*lam1H + 256*pow(pi,2)*lam1H - 8*Lb*(24*pow(lam1H,2) + 2*pow(lam3H,2) + 2*lam3H*lam4H + pow(lam4H,2) + pow(lam5H,2))))/(256.*pow(pi,2));
      std::complex<double> lam2H3d = -0.00390625*(T*(pow(g1,4)*(-2 + 3*Lb) + pow(g2,4)*(-6 + 9*Lb) - 72*pow(g2,2)*Lb*lam2H + pow(g1,2)*(pow(g2,2)*(-4 + 6*Lb) - 24*Lb*lam2H) + 8*(-32*pow(pi,2)*lam2H - 6*Lf*(pow(yt2,4) - 2*pow(yt2,2)*lam2H) + Lb*(24*pow(lam2H,2) + 2*pow(lam3H,2) + 2*lam3H*lam4H + pow(lam4H,2) + pow(lam5H,2)))))/pow(pi,2);
      std::complex<double> lam3H3d = (T*((pow(g1,4) - 2*pow(g1,2)*pow(g2,2) + 3*pow(g2,4))*(2 - 3*Lb) + 128*pow(pi,2)*lam3H + 12*(pow(g1,2)*Lb + 3*pow(g2,2)*Lb - 2*Lf*pow(yt2,2))*lam3H - 8*Lb*(2*lam3H*(3*(lam1H + lam2H) + lam3H) + 2*(lam1H + lam2H)*lam4H + pow(lam4H,2) + pow(lam5H,2))))/(128.*pow(pi,2));
      std::complex<double> lam4H3d = -0.03125*(T*(-9*pow(g2,2)*Lb*lam4H - 32*pow(pi,2)*lam4H + 6*Lf*pow(yt2,2)*lam4H + 4*Lb*lam1H*lam4H + 4*Lb*lam2H*lam4H + 8*Lb*lam3H*lam4H + 4*Lb*pow(lam4H,2) + pow(g1,2)*(pow(g2,2)*(-2 + 3*Lb) - 3*Lb*lam4H) + 8*Lb*pow(lam5H,2)))/pow(pi,2);
      std::complex<double> lam5H3d = (T*(3*pow(g1,2)*Lb + 9*pow(g2,2)*Lb + 32*pow(pi,2) - 6*Lf*pow(yt2,2) - 4*Lb*(lam1H + lam2H + 2*lam3H + 3*lam4H))*lam5H)/(32.*pow(pi,2));
      std::complex<double> lam6H3d = 0;
      std::complex<double> lam7H3d = 0;
      //---------------------------------------------------------------------------------
      std::complex<double> lamVLL1 = (-181*pow(g1,4)*T)/(36.*pow(pi,2));
      std::complex<double> lamVLL2 = -0.25*(pow(g1,2)*pow(g2,2)*T)/pow(pi,2);
      std::complex<double> lamVLL3 = (pow(g2,4)*T)/(4.*pow(pi,2));
      std::complex<double> lamVLL4 = (-11*pow(g1,2)*pow(g3,2)*T)/(12.*pow(pi,2));
      std::complex<double> lamVLL5 = (-3*pow(g2,2)*pow(g3,2)*T)/(4.*pow(pi,2));
      std::complex<double> lamVLL6 = -0.25*(g1*pow(g3,3)*T)/pow(pi,2);
      std::complex<double> lamVLL7 = (pow(g3,4)*T)/(4.*pow(pi,2));
      std::complex<double> lamVL1 = -0.25*(pow(g3,2)*T*pow(yt2,2))/pow(pi,2);
      std::complex<double> lamVL2 = -0.005208333333333333*(g1*g2*T*(pow(g2,2)*(-5 - 21*Lb + 12*Lf) + pow(g1,2)*(-21 + Lb + 20*Lf) - 12*(8*pow(pi,2) + 2*lam1H + lam4H)))/pow(pi,2);
      std::complex<double> lamVL3 = -0.005208333333333333*(g1*g2*T*(pow(g2,2)*(-5 - 21*Lb + 12*Lf) + pow(g1,2)*(-21 + Lb + 20*Lf) - 12*(8*pow(pi,2) + pow(yt2,2) + 2*lam2H + lam4H)))/pow(pi,2);
      std::complex<double> lamVL4 = -0.005208333333333333*(pow(g1,2)*T*(-9*pow(g2,2) + pow(g1,2)*(-39 + 2*Lb + 40*Lf) - 12*(8*pow(pi,2) + 6*lam1H + 2*lam3H + lam4H)))/pow(pi,2);
      std::complex<double> lamVL5 = (pow(g2,2)*T*(3*pow(g1,2) + pow(g2,2)*(73 + 42*Lb - 24*Lf) + 12*(8*pow(pi,2) + 6*lam1H + 2*lam3H + lam4H)))/(192.*pow(pi,2));
      std::complex<double> lamVL6 = -0.005208333333333333*(pow(g1,2)*T*(-9*pow(g2,2) + pow(g1,2)*(-39 + 2*Lb + 40*Lf) - 96*pow(pi,2) + 68*pow(yt2,2) - 12*(6*lam2H + 2*lam3H + lam4H)))/pow(pi,2);
      std::complex<double> lamVL7 = (pow(g2,2)*T*(3*pow(g1,2) + pow(g2,2)*(73 + 42*Lb - 24*Lf) + 12*(8*pow(pi,2) - 3*pow(yt2,2) + 6*lam2H + 2*lam3H + lam4H)))/(192.*pow(pi,2));
      std::complex<double> lamVL8 = 0;
      std::complex<double> lamVL9 = 0;
      std::complex<double> lamVL10 = 0;
      //---------------------------------------------------------------------------------
      std::complex<double> musqSU2 = 2*pow(g2,2)*pow(T,2) + (pow(g2,2)*(72*m1 + 72*m2 + \
        3*pow(T,2)*(-3*pow(g1,2) - 72*pow(g3,2) - \
        3*pow(yt2,2) + 4*(3*lam1H + 3*lam2H + 2*lam3H + lam4H) + \
        pow(g2,2)*(115 + 144*EulerGamma - 672*log(2))) + \
        8*pow(g2,2)*pow(T,2)*(-54*log(pi) - 25*log(T) + \
        25*log(mu4) + 29*log(mu4/T))))/(576.*pow(pi,2));

      std::complex<double> musqSU3 = (pow(g3,2)*pow(T,2)*(384 + (-11*pow(g1,2) + \
              3*(-9*pow(g2,2) - 4*pow(yt2,2) + 8*pow(g3,2)*(5 + \
              14*EulerGamma - 22*log(2))) + \
              24*pow(g3,2)*(-3*log(pi*T) + 3*log(mu4) + \
              11*log(mu4/(4.*pi*T))))/pow(pi,2)))/192.;

      std::complex<double> musqU1 = 2*pow(g1,2)*pow(T,2) + (pow(g1,2)*(72*m1 + 72*m2 + \
              pow(T,2)*(3*(-9*pow(g2,2) - 88*pow(g3,2) - \
              11*pow(yt2,2) + 4*(3*lam1H + 3*lam2H + 2*lam3H + lam4H)) + \
              pow(g1,2)*(251 - 1008*EulerGamma + \
              48*log(4*pow(pi,21)))) + \
              8*pow(g1,2)*pow(T,2)*(101*(log(T) - log(mu4)) - \
              25*log(mu4/T))))/(576.*pow(pi,2));
      //---------------------------------------------------------------------------------
      mu3 = T;
      //mu3dUS = sqrt(g1sq) * scaleFactor3;
      //---------------------------------------------------------------------------------
      std::complex<double> m12R3d = (m12R*(3*pow(g1,2)*Lb + 9*pow(g2,2)*Lb - \
              2*(-32*pow(pi,2) + 3*Lf*pow(yt2,2) + 2*Lb*lam3H + \
              4*Lb*lam4H + 6*Lb*lam5H)) + 2*log(mu3/mu4)*(-12*lam1H3d*lam6H3d - \
              6*lam3H3d*lam6H3d - 6*lam4H3d*lam6H3d - 6*lam5H3d*lam6H3d - \
              12*lam2H3d*lam7H3d - 6*lam3H3d*lam7H3d - 6*lam4H3d*lam7H3d - \
              6*lam5H3d*lam7H3d + 3*g13d2*(lam6H3d + lam7H3d) - lamVL4*lamVL8 - \
              lamVL6*lamVL8 + 6*lamVL2*lamVL9 + 6*lamVL3*lamVL9 - \
              3*lamVL5*lamVL10 - 3*lamVL7*lamVL10 + 3*g23d2*(3*lam6H3d + \
              3*lam7H3d + 4*lamVL10)))/(64.*pow(pi,2));
      
      std::complex<double> m13d = m1 + (pow(T,2)*(3*pow(g1,2) + 9*pow(g2,2) + 4*(6*lam1H \
              + 2*lam3H + lam4H)))/48. + (-96*Lb*(6*m1*lam1H + m2*(2*lam3H + \
              lam4H)) + pow(g1,2)*(12*Lb*(6*m1 + \
              pow(T,2)*(6*pow(g2,2) - 6*lam1H - 2*lam3H - lam4H)) + \
              2*pow(T,2)*(2*(6*lam1H + 2*lam3H + lam4H)*(1 + 6*EulerGamma - \
              72*log(Glaisher)) - 9*pow(g2,2)*(1 + 5*EulerGamma - \
              60*log(Glaisher)))) + 6*pow(g2,2)*(36*Lb*m1 + \
              3*Lb*pow(T,2)*(6*lam1H + 2*lam3H + lam4H) + \
              pow(T,2)*(6*lam1H + 2*lam3H + lam4H)*(2 + 12*EulerGamma - 9*Lb - \
              144*log(Glaisher))) + 3*pow(g2,4)*pow(T,2)*(67 + \
              75*EulerGamma - 84*Lb + 12*Lf - 900*log(Glaisher)) - \
              pow(g1,4)*pow(T,2)*(-17 + 27*EulerGamma + 50*Lb - 20*Lf - \
              324*log(Glaisher)) - \
              4*pow(T,2)*(9*Lb*pow(yt2,2)*(2*lam3H + lam4H) - \
              3*Lf*pow(yt2,2)*(2*lam3H + lam4H) + 2*Lb*(-2*pow(lam3H,2) - \
              2*lam3H*lam4H - 5*pow(lam4H,2) + 6*lam1H*(2*lam3H + lam4H) + \
              6*lam2H*(2*lam3H + lam4H) - 9*pow(lam5H,2)) + \
              12*EulerGamma*(12*pow(lam1H,2) + 2*pow(lam3H,2) + \
              2*lam3H*lam4H + 2*pow(lam4H,2) + 3*pow(lam5H,2)) - \
              144*(12*pow(lam1H,2) + 2*pow(lam3H,2) + 2*lam3H*lam4H + \
              2*pow(lam4H,2) + 3*pow(lam5H,2))*log(Glaisher)) + \
              6*log(mu3/mu4)*(7*(pow(g13d2,2.)) - 33*(pow(g23d2,2.)) + 2*g13d2*(9*g23d2 - \
              4*(6*lam1H3d + 2*lam3H3d + lam4H3d)) - 24*g23d2*(6*lam1H3d + \
              2*lam3H3d + lam4H3d + 4*lamVL5) + 8*(24*pow(lam1H3d,2) + \
              4*pow(lam3H3d,2) + 4*lam3H3d*lam4H3d + 4*pow(lam4H3d,2) + \
              6*pow(lam5H3d,2) + 18*pow(lam6H3d,2) + \
              6*pow(lam7H3d,2) + 6*pow(lamVL2,2) + pow(lamVL4,2) \
              + 3*pow(lamVL5,2) + pow(lamVL8,2) + \
              6*pow(lamVL9,2) + \
              3*pow(lamVL10,2))))/(1536.*pow(pi,2));

      std::complex<double> m23d = (9*pow(g2,2)*(72*Lb*m2 + 3*Lb*pow(T,2)*(7*pow(yt2,2) - \
              4*(6*lam2H + 2*lam3H + lam4H)) + pow(T,2)*(96*pow(pi,2) \
              - 3*(2 + Lf)*pow(yt2,2) + 4*(6*lam2H + 2*lam3H + lam4H)*(1 + \
              6*EulerGamma - 72*log(Glaisher)))) + \
              9*pow(g2,4)*pow(T,2)*(67 + 75*EulerGamma - 84*Lb + 12*Lf - \
              900*log(Glaisher)) + 3*pow(g1,4)*pow(T,2)*(17 - \
              27*EulerGamma - 50*Lb + 20*Lf + 324*log(Glaisher)) - \
              12*(16*pow(g3,2)*(3 + Lb - 4*Lf)*pow(T,2)*pow(yt2,2) - \
              96*pow(pi,2)*pow(T,2)*pow(yt2,2) - \
              9*Lb*pow(T,2)*pow(yt2,4) - \
              192*pow(pi,2)*pow(T,2)*lam2H + \
              54*Lb*pow(T,2)*pow(yt2,2)*lam2H + \
              18*Lf*pow(T,2)*pow(yt2,2)*lam2H + \
              144*EulerGamma*pow(T,2)*pow(lam2H,2) - \
              24*m2*(16*pow(pi,2) - 3*Lf*pow(yt2,2) - 6*Lb*lam2H) + \
              48*Lb*m1*lam3H - 64*pow(pi,2)*pow(T,2)*lam3H + \
              12*Lf*pow(T,2)*pow(yt2,2)*lam3H + \
              24*Lb*pow(T,2)*lam1H*lam3H + 24*Lb*pow(T,2)*lam2H*lam3H + \
              24*EulerGamma*pow(T,2)*pow(lam3H,2) - \
              4*Lb*pow(T,2)*pow(lam3H,2) + 24*Lb*m1*lam4H - \
              32*pow(pi,2)*pow(T,2)*lam4H + \
              6*Lf*pow(T,2)*pow(yt2,2)*lam4H + \
              12*Lb*pow(T,2)*lam1H*lam4H + 12*Lb*pow(T,2)*lam2H*lam4H + \
              24*EulerGamma*pow(T,2)*lam3H*lam4H - \
              4*Lb*pow(T,2)*lam3H*lam4H + \
              24*EulerGamma*pow(T,2)*pow(lam4H,2) - \
              10*Lb*pow(T,2)*pow(lam4H,2) + \
              36*EulerGamma*pow(T,2)*pow(lam5H,2) - \
              18*Lb*pow(T,2)*pow(lam5H,2) - \
              1728*pow(T,2)*pow(lam2H,2)*log(Glaisher) - \
              288*pow(T,2)*pow(lam3H,2)*log(Glaisher) - \
              288*pow(T,2)*lam3H*lam4H*log(Glaisher) - \
              288*pow(T,2)*pow(lam4H,2)*log(Glaisher) - \
              432*pow(T,2)*pow(lam5H,2)*log(Glaisher)) + \
              pow(g1,2)*(Lb*(216*m2 + pow(T,2)*(216*pow(g2,2) + \
              47*pow(yt2,2) - 36*(6*lam2H + 2*lam3H + lam4H))) + \
              pow(T,2)*(288*pow(pi,2) - 66*pow(yt2,2) + \
              55*Lf*pow(yt2,2) + 72*lam2H + 432*EulerGamma*lam2H + 24*lam3H + \
              144*EulerGamma*lam3H + 12*lam4H + 72*EulerGamma*lam4H - \
              54*pow(g2,2)*(1 + 5*EulerGamma - 60*log(Glaisher)) - \
              5184*lam2H*log(Glaisher) - 1728*lam3H*log(Glaisher) - \
              864*lam4H*log(Glaisher))) + 18*log(mu3/mu4)*(7*(pow(g13d2,2.)) - \
              33*(pow(g23d2,2.)) + 2*g13d2*(9*g23d2 - 4*(6*lam2H3d + 2*lam3H3d + \
              lam4H3d)) - 24*g23d2*(6*lam2H3d + 2*lam3H3d + lam4H3d + 4*lamVL7) + \
              8*(24*pow(lam2H3d,2) + 4*pow(lam3H3d,2) + 4*lam3H3d*lam4H3d \
              + 4*pow(lam4H3d,2) + 6*pow(lam5H3d,2) + \
              6*pow(lam6H3d,2) + 18*pow(lam7H3d,2) - 48*g33d2*lamVL1 + \
              8*pow(lamVL1,2) + 6*pow(lamVL3,2) + \
              pow(lamVL6,2) + 3*pow(lamVL7,2) + pow(lamVL8,2) \
              + 6*pow(lamVL9,2) + \
              3*pow(lamVL10,2))))/(4608.*pow(pi,2));
      //---------------------------------------------------------------------------------
      /*if ( matching_flag == 0 ) {
        return {g1sq*T, g2sq*T, g3sq*T, b1/sqrt(T), mphiSq3d_LO, msSq3d_LO, a1*sqrt(T), b3*sqrt(T), lamH*T, lamHS*T, lamS*T};
      } else if ( matching_flag == 1 ) {
        return {g1sq3d, g2sq3d, g3sq3d, b13d_NLO, mphiSq3d, msSq3d, a13d, b33d, lamH3d, lamHS3d, lamS3d};
      }*/
      //---------------------------------------------------------------------------------
      std::complex<double> lam1H3dUS = lam1H3d - ((4*pow(lamVL2,2))/(sqrt(musqSU2) + sqrt(musqU1)) + pow(lamVL4,2)/sqrt(musqU1) + (3*pow(lamVL5,2))/sqrt(musqSU2))/(32.*pi);
      std::complex<double> lam2H3dUS = lam2H3d - ((8*pow(lamVL1,2))/sqrt(musqSU3) + (4*pow(lamVL3,2))/(sqrt(musqSU2) + sqrt(musqU1)) + pow(lamVL6,2)/sqrt(musqU1) + (3*pow(lamVL7,2))/sqrt(musqSU2))/(32.*pi);
      std::complex<double> lam3H3dUS = (16*pi*lam3H3d*(musqSU2*sqrt(musqU1) + sqrt(musqSU2)*musqU1) - musqSU2*lamVL4*lamVL6 - 3*musqU1*lamVL5*lamVL7 + sqrt(musqSU2)*sqrt(musqU1)*(4*lamVL2*lamVL3 - lamVL4*lamVL6 - 3*lamVL5*lamVL7 - 8*pow(lamVL9,2)))/(16.*pi*(musqSU2*sqrt(musqU1) + sqrt(musqSU2)*musqU1));
      std::complex<double> lam4H3dUS = -0.0625*(-16*pi*lam4H3d*(musqSU2*sqrt(musqU1) + sqrt(musqSU2)*musqU1) + musqSU2*pow(lamVL8,2) + 3*musqU1*pow(lamVL10,2) + sqrt(musqSU2)*sqrt(musqU1)*(8*lamVL2*lamVL3 + pow(lamVL8,2) - 4*pow(lamVL9,2) + 3*pow(lamVL10,2)))/(pi*(musqSU2*sqrt(musqU1) + sqrt(musqSU2)*musqU1));
      std::complex<double> lam5H3dUS = lam5H3d - (pow(lamVL8,2)/sqrt(musqU1) + (4*pow(lamVL9,2))/(sqrt(musqSU2) + sqrt(musqU1)) + (3*pow(lamVL10,2))/sqrt(musqSU2))/(16.*pi);
      std::complex<double> lam6H3dUS = lam6H3d - ((lamVL4*lamVL8)/sqrt(musqU1) - (4*lamVL2*lamVL9)/(sqrt(musqSU2) + sqrt(musqU1)) + (3*lamVL5*lamVL10)/sqrt(musqSU2))/(16.*pi);
      std::complex<double> lam7H3dUS = lam7H3d - ((lamVL6*lamVL8)/sqrt(musqU1) - (4*lamVL3*lamVL9)/(sqrt(musqSU2) + sqrt(musqU1)) + (3*lamVL7*lamVL10)/sqrt(musqSU2))/(16.*pi);
      std::complex<double> g13dUS2 = g13d2;
      std::complex<double> g23dUS2 = g23d2 - pow(g23d2,2.)/(24.*pi*sqrt(musqSU2));
      std::complex<double> g33dUS2 = g33d2 - pow(g33d2,2.)/(16.*pi*sqrt(musqSU3));
      //---------------------------------------------------------------------------------
      std::complex<double> m12R3dUS = m12R3d + (sqrt(musqU1)*lamVL8 + \
              3*sqrt(musqSU2)*lamVL10)/(8.*pi) - (-2*(1 + \
              2*log(mu3/(2.*sqrt(musqU1))))*(lamVL4 + lamVL6)*lamVL8 + \
              12*(1 + 2*log(mu3/(sqrt(musqSU2) + sqrt(musqU1))))*(lamVL2 \
              + lamVL3)*lamVL9 + 12*g23d2*(1 + \
              4*log(mu3/(2.*sqrt(musqSU2))))*lamVL10 - 6*(1 + \
              2*log(mu3/(2.*sqrt(musqSU2))))*(lamVL5 + lamVL7)*lamVL10 \
              + (lamVL8*(sqrt(musqU1)*lamVLL1 + 3*sqrt(musqSU2)*lamVLL2 \
              + 8*sqrt(musqSU3)*lamVLL4))/sqrt(musqU1) + \
              (3*lamVL10*(sqrt(musqU1)*lamVLL2 + \
              5*sqrt(musqSU2)*lamVLL3 + \
              8*sqrt(musqSU3)*lamVLL5))/sqrt(musqSU2))/(128.*pow(pi,\
              2));
      std::complex<double> m13dUS = m13d - (sqrt(musqU1)*lamVL4 + \
              3*sqrt(musqSU2)*lamVL5)/(8.*pi) + \
              (-6*pow(g23d2,2.)*log(mu3/(2.*sqrt(musqSU2))) + 12*g23d2*(1 + \
              4*log(mu3/(2.*sqrt(musqSU2))))*lamVL5 - 2*(1 + \
              2*log(mu3/(2.*sqrt(musqU1))))*(pow(lamVL4,2) + \
              pow(lamVL8,2)) - 12*(1 + 2*log(mu3/(sqrt(musqSU2) + \
              sqrt(musqU1))))*(pow(lamVL2,2) + pow(lamVL9,2)) - \
              6*(1 + 2*log(mu3/(2.*sqrt(musqSU2))))*(pow(lamVL5,2) + \
              pow(lamVL10,2)) + (lamVL4*(sqrt(musqU1)*lamVLL1 + \
              3*sqrt(musqSU2)*lamVLL2 + \
              8*sqrt(musqSU3)*lamVLL4))/sqrt(musqU1) + \
              (3*lamVL5*(sqrt(musqU1)*lamVLL2 + 5*sqrt(musqSU2)*lamVLL3 \
              + 8*sqrt(musqSU3)*lamVLL5))/sqrt(musqSU2))/(128.*pow(\
              pi,2));
      std::complex<double> m23dUS = m23d - (8*sqrt(musqSU3)*lamVL1 + sqrt(musqU1)*lamVL6 + \
              3*sqrt(musqSU2)*lamVL7)/(8.*pi) + \
              (-6*pow(g23d2,2.)*log(mu3/(2.*sqrt(musqSU2))) + 48*g33d2*(1 + \
              4*log(mu3/(2.*sqrt(musqSU3))))*lamVL1 - 16*(1 + \
              2*log(mu3/(2.*sqrt(musqSU3))))*pow(lamVL1,2) + \
              12*g23d2*(1 + 4*log(mu3/(2.*sqrt(musqSU2))))*lamVL7 - 2*(1 + \
              2*log(mu3/(2.*sqrt(musqU1))))*(pow(lamVL6,2) + \
              pow(lamVL8,2)) - 12*(1 + 2*log(mu3/(sqrt(musqSU2) + \
              sqrt(musqU1))))*(pow(lamVL3,2) + pow(lamVL9,2)) - \
              6*(1 + 2*log(mu3/(2.*sqrt(musqSU2))))*(pow(lamVL7,2) + \
              pow(lamVL10,2)) + (lamVL6*(sqrt(musqU1)*lamVLL1 + \
              3*sqrt(musqSU2)*lamVLL2 + \
              8*sqrt(musqSU3)*lamVLL4))/sqrt(musqU1) + \
              (3*lamVL7*(sqrt(musqU1)*lamVLL2 + 5*sqrt(musqSU2)*lamVLL3 \
              + 8*sqrt(musqSU3)*lamVLL5))/sqrt(musqSU2) + \
              (8*lamVL1*(sqrt(musqU1)*lamVLL4 + 3*sqrt(musqSU2)*lamVLL5 \
              + 10*sqrt(musqSU3)*lamVLL7))/sqrt(musqSU3))/(128.*pow(\
              pi,2));
      //---------------------------------------------------------------------------------
      std::complex<double> mu3US = mu3;
      //---------------------------------------------------------------------------------
      
      return {g13dUS2, g23dUS2, g33dUS2, lam1H3dUS, lam2H3dUS, lam3H3dUS, lam4H3dUS, lam5H3dUS, lam6H3dUS, lam7H3dUS, m12R3dUS, m13dUS, m23dUS, mu4, mu3US};
      
    }

  private :
    double pi = 3.141459;
    std::vector<double> LagrangianParams;
    std::vector<double> InputParams;

    /*const double v = SM::v;
    const double Mh = SM::mh; // Captial for physical mass.
    const double MhSq = Mh*Mh;
    const double g1 = SM::gp;
    const double g2 = SM::g;*/
  
    // const double scaleFactor = M_PI;

};

} // namespace EffectivePotential

#endif
