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

#include <vector>
#include <cmath>
#include <complex>
#include "ComplexOperators.hpp"

class DR_2HDM_RenormParam {

  private:
    
    double real(std::complex<double> comp){
        return comp.real();
    }

    double pi = 3.14159;

    double alphaEM = 1./137.035999679;

    std::complex<double> alphaS = 0.1181;
    std::complex<double> gs = sqrt(4*pi*alphaS);

    double Nf = 3;

    double mt = 173.1;
    double mc = 1.28;
    double mup = 2.2*pow(10,-3);
    double mb = 4.18;
    double ms = 96*pow(10,-3);
    double md = 4.7*pow(10,-3);

    double mtau = 1.77686;
    double mmu = 105.6583745*pow(10,-3);
    double me = 0.5109989461*pow(10,-3);
    double mNuTau = 0;
    double mNuMu = 0;
    double mNuE = 0;

    double vev = 246;
    double mW = 80.385;
    double mZ = 91.1876;

    std::complex<double> Mt = mt;
    std::complex<double> Mc = mc;
    std::complex<double> Mu = mup;
    std::complex<double> Mb = mb;
    std::complex<double> Ms = ms;
    std::complex<double> Md = md;

    std::complex<double> Mtau = mtau;
    std::complex<double> Mmu = mmu;
    std::complex<double> Me = me;
    std::complex<double> MNuTau = mNuTau+ pow(10.,-5);
    std::complex<double> MNuMu = mNuMu+ pow(10.,-5);
    std::complex<double> MNuE = mNuE + pow(10.,-5);

    std::complex<double> MW = mW;
    std::complex<double> MZ = mZ;

    std::complex<double> Mh, MH, MA, MHpm, mu, alpha, beta, rgLam;

  public:

    DR_2HDM_RenormParam(std::vector<double> params) {
        Mh = params[0];
        MH = params[1];
        MA = params[2] + pow(10., -6.);
        MHpm = params[3];
        mu = params[4];
        alpha = params[5];
        beta = params[6];
        rgLam = params[7];
    }
    
    std::vector<double> OneLoopRenormalised(){
        
        std::vector<double> Params4D = {0.130, 0.418, 1.489, 0.3, 0.0925, 3.675, -1.78, -1.792, sqrt(0.998), -6400, 12942, -2751};

        return Params4D;

    }


  };