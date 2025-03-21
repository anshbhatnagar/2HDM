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

class DR_2HDM_Potential {

  public:
    const double pi = 3.14159;

    std::complex<double> V0( Eigen::VectorXd phi, double T, std::vector<std::complex<double>> par ) {

      std::complex<double> g13dUS = sqrt(par[0]);
      std::complex<double> g23dUS = sqrt(par[1]);
      std::complex<double> g33dUS = sqrt(par[2]);
      std::complex<double> lam1H3dUS = par[3];
      std::complex<double> lam2H3dUS = par[4];
      std::complex<double> lam3H3dUS = par[5];
      std::complex<double> lam4H3dUS = par[6];
      std::complex<double> lam5H3dUS = par[7];
      std::complex<double> lam6H3dUS = par[8];
      std::complex<double> lam7H3dUS = par[9];
      std::complex<double> m12R3dUS = par[10];
      std::complex<double> m13dUS = par[11];
      std::complex<double> m23dUS = par[12];
      std::complex<double> mu = par[13];
      std::complex<double> mu3US = par[14];
    
      std::complex<double> phi1(phi[0]/sqrt(T + 1e-15), 0);
      std::complex<double> phi2(phi[1]/sqrt(T + 1e-15), 0);

      return (phi1*(m13dUS*phi1 - m12R3dUS*phi2) + phi2*(-(m12R3dUS*phi1) + \
        m23dUS*phi2))/2. + (phi2*(2*(lam3H3dUS + lam4H3dUS + \
        lam5H3dUS)*pow(phi1,2)*phi2 + phi2*((lam3H3dUS + lam4H3dUS + \
        lam5H3dUS)*pow(phi1,2) + 6*lam2H3dUS*pow(phi2,2))) + \
        phi1*(2*(lam3H3dUS + lam4H3dUS + lam5H3dUS)*phi1*pow(phi2,2) + \
        phi1*(6*lam1H3dUS*pow(phi1,2) + (lam3H3dUS + lam4H3dUS + \
        lam5H3dUS)*pow(phi2,2))))/24.;
    }

    std::complex<double> V1( Eigen::VectorXd phi, double T, std::vector<std::complex<double>> par ) {

      std::complex<double> g13dUS = sqrt(par[0]);
      std::complex<double> g23dUS = sqrt(par[1]);
      std::complex<double> g33dUS = sqrt(par[2]);
      std::complex<double> lam1H3dUS = par[3];
      std::complex<double> lam2H3dUS = par[4];
      std::complex<double> lam3H3dUS = par[5];
      std::complex<double> lam4H3dUS = par[6];
      std::complex<double> lam5H3dUS = par[7];
      std::complex<double> lam6H3dUS = par[8];
      std::complex<double> lam7H3dUS = par[9];
      std::complex<double> m12R3dUS = par[10];
      std::complex<double> m13dUS = par[11];
      std::complex<double> m23dUS = par[12];
      std::complex<double> mu = par[13];
      std::complex<double> mu3US = par[14];
    
      std::complex<double> phi1(phi[0]/sqrt(T + 1e-15), 0);
      std::complex<double> phi2(phi[1]/sqrt(T + 1e-15), 0);

      std::complex<double> Compile_2 = 1./pi;
      std::complex<double> Compile_3 = pow(g23dUS,2);
      std::complex<double> Compile_5 = pow(phi1,2);
      std::complex<double> Compile_6 = pow(phi2,2);
      std::complex<double> Compile_7 = Compile_5 + Compile_6;
      std::complex<double> Compile_21 = 2*m13dUS;
      std::complex<double> Compile_22 = -2*m23dUS;
      std::complex<double> Compile_30 = -lam3H3dUS;
      std::complex<double> Compile_29 = 2*lam1H3dUS;
      std::complex<double> Compile_31 = Compile_29 + Compile_30;
      std::complex<double> Compile_35 = lam4H3dUS + lam5H3dUS;
      std::complex<double> Compile_42 = -2*lam2H3dUS;
      std::complex<double> Compile_43 = lam3H3dUS + Compile_42;
      std::complex<double> Compile_23 = 2*lam1H3dUS*Compile_5;
      std::complex<double> Compile_26 = lam3H3dUS*Compile_6;
      std::complex<double> Compile_24 = -(lam3H3dUS*Compile_5);
      std::complex<double> Compile_25 = -2*lam2H3dUS*Compile_6;
      std::complex<double> Compile_27 = pow(m12R3dUS,2);
      std::complex<double> Compile_28 = 16*Compile_27;
      std::complex<double> Compile_32 = Compile_31*Compile_5;
      std::complex<double> Compile_33 = Compile_21 + Compile_22 + Compile_32;
      std::complex<double> Compile_34 = pow(Compile_33,2);
      std::complex<double> Compile_36 = -16*m12R3dUS*phi1*phi2*Compile_35;
      std::complex<double> Compile_37 = -m23dUS;
      std::complex<double> Compile_38 = m13dUS + Compile_37;
      std::complex<double> Compile_39 = 2*lam2H3dUS;
      std::complex<double> Compile_40 = Compile_30 + Compile_39;
      std::complex<double> Compile_41 = -2*Compile_38*Compile_40;
      std::complex<double> Compile_44 = Compile_31*Compile_43;
      std::complex<double> Compile_45 = pow(Compile_35,2);
      std::complex<double> Compile_46 = 2*Compile_45;
      std::complex<double> Compile_47 = Compile_44 + Compile_46;
      std::complex<double> Compile_48 = Compile_47*Compile_5;
      std::complex<double> Compile_49 = Compile_41 + Compile_48;
      std::complex<double> Compile_50 = 2*Compile_49*Compile_6;
      std::complex<double> Compile_51 = pow(Compile_43,2);
      std::complex<double> Compile_52 = pow(phi2,4);
      std::complex<double> Compile_53 = Compile_51*Compile_52;
      std::complex<double> Compile_54 = Compile_28 + Compile_34 + \
      Compile_36 + Compile_50 + Compile_53;
      std::complex<double> Compile_55 = sqrt(Compile_54);
      std::complex<double> Compile_56 = Compile_21 + Compile_22 + \
      Compile_23 + Compile_24 + Compile_25 + Compile_26 + Compile_55;
      std::complex<double> Compile_18 = 4*m23dUS;
      std::complex<double> Compile_19 = 2*lam3H3dUS*Compile_5;
      std::complex<double> Compile_20 = 4*lam2H3dUS*Compile_6;
      std::complex<double> Compile_58 = -2*m12R3dUS;
      std::complex<double> Compile_59 = phi1*phi2*Compile_35;
      std::complex<double> Compile_60 = Compile_58 + Compile_59;
      std::complex<double> Compile_61 = pow(Compile_60,-2);
      std::complex<double> Compile_79 = -2*m13dUS;
      std::complex<double> Compile_80 = 2*m23dUS;
      std::complex<double> Compile_81 = -2*lam1H3dUS*Compile_5;
      std::complex<double> Compile_83 = 2*lam2H3dUS*Compile_6;
      std::complex<double> Compile_66 = 2*m12R3dUS;
      std::complex<double> Compile_67 = -(phi1*phi2*Compile_35);
      std::complex<double> Compile_68 = Compile_66 + Compile_67;
      std::complex<double> Compile_69 = 1./Compile_68;
      std::complex<double> Compile_90 = -phi2;
      std::complex<double> Compile_91 = phi1 + Compile_90;
      std::complex<double> Compile_92 = phi1 + phi2;
      std::complex<double> Compile_93 = lam3H3dUS*Compile_91*Compile_92;
      std::complex<double> Compile_94 = Compile_55 + Compile_79 + \
      Compile_80 + Compile_81 + Compile_83 + Compile_93;
      std::complex<double> Compile_106 = -lam5H3dUS;
      std::complex<double> Compile_107 = lam3H3dUS + lam4H3dUS + Compile_106;
      std::complex<double> Compile_118 = -lam4H3dUS;
      std::complex<double> Compile_111 = -(lam5H3dUS*phi1*phi2);
      std::complex<double> Compile_112 = m12R3dUS + Compile_111;
      std::complex<double> Compile_87 = lam1H3dUS*Compile_5;
      std::complex<double> Compile_117 = Compile_107*Compile_91*Compile_92;
      std::complex<double> Compile_119 = lam5H3dUS + Compile_118 + \
      Compile_29 + Compile_30;
      std::complex<double> Compile_120 = Compile_119*Compile_5;
      std::complex<double> Compile_121 = Compile_120 + Compile_21 + \
      Compile_22;
      std::complex<double> Compile_122 = pow(Compile_121,2);
      std::complex<double> Compile_123 = -32*m12R3dUS*lam5H3dUS*phi1*phi2;
      std::complex<double> Compile_124 = lam5H3dUS + Compile_118 + \
      Compile_30 + Compile_39;
      std::complex<double> Compile_125 = 2*Compile_124*Compile_38;
      std::complex<double> Compile_126 = Compile_118 + Compile_29 + \
      Compile_30;
      std::complex<double> Compile_127 = lam3H3dUS + lam4H3dUS + \
      Compile_42;
      std::complex<double> Compile_128 = -(Compile_126*Compile_127);
      std::complex<double> Compile_129 = lam1H3dUS + lam2H3dUS + \
      Compile_118 + Compile_30;
      std::complex<double> Compile_130 = 2*lam5H3dUS*Compile_129;
      std::complex<double> Compile_131 = pow(lam5H3dUS,2);
      std::complex<double> Compile_132 = -7*Compile_131;
      std::complex<double> Compile_133 = Compile_128 + Compile_130 + \
      Compile_132;
      std::complex<double> Compile_134 = Compile_133*Compile_5;
      std::complex<double> Compile_135 = Compile_125 + Compile_134;
      std::complex<double> Compile_136 = -2*Compile_135*Compile_6;
      std::complex<double> Compile_137 = lam3H3dUS + lam4H3dUS + \
      Compile_106 + Compile_42;
      std::complex<double> Compile_138 = pow(Compile_137,2);
      std::complex<double> Compile_139 = Compile_138*Compile_52;
      std::complex<double> Compile_140 = Compile_122 + Compile_123 + \
      Compile_136 + Compile_139 + Compile_28;
      std::complex<double> Compile_141 = sqrt(Compile_140);
      std::complex<double> Compile_142 = Compile_117 + Compile_141 + \
      Compile_79 + Compile_80 + Compile_81 + Compile_83;
      std::complex<double> Compile_113 = 1./Compile_112;
      std::complex<double> Compile_108 = Compile_107*Compile_5;
      std::complex<double> Compile_109 = Compile_108 + Compile_80 + \
      Compile_83;
      std::complex<double> Compile_110 = 8*Compile_109;
      std::complex<double> Compile_144 = pow(Compile_112,-2);
      std::complex<double> Compile_145 = (Compile_107*Compile_6)/2.;
      std::complex<double> Compile_146 = m13dUS + Compile_145 + Compile_87;
      std::complex<double> Compile_158 = \
      -(Compile_107*Compile_91*Compile_92);
      std::complex<double> Compile_159 = Compile_141 + Compile_158 + \
      Compile_21 + Compile_22 + Compile_23 + Compile_25;
      std::complex<double> Compile_82 = lam3H3dUS*Compile_5;
      std::complex<double> Compile_178 = 6*lam2H3dUS*Compile_6;
      std::complex<double> Compile_84 = -(lam3H3dUS*Compile_6);
      std::complex<double> Compile_176 = lam3H3dUS + lam4H3dUS + lam5H3dUS;
      std::complex<double> Compile_192 = 6*lam2H3dUS;
      std::complex<double> Compile_195 = -6*lam2H3dUS;
      std::complex<double> Compile_196 = lam3H3dUS + lam4H3dUS + lam5H3dUS \
      + Compile_195;
      std::complex<double> Compile_181 = -6*lam1H3dUS*Compile_5;
      std::complex<double> Compile_186 = 6*lam1H3dUS;
      std::complex<double> Compile_187 = Compile_106 + Compile_118 + \
      Compile_186 + Compile_30;
      std::complex<double> Compile_188 = Compile_187*Compile_5;
      std::complex<double> Compile_189 = Compile_188 + Compile_21 + \
      Compile_22;
      std::complex<double> Compile_190 = pow(Compile_189,2);
      std::complex<double> Compile_191 = \
      -32*m12R3dUS*phi1*phi2*Compile_176;
      std::complex<double> Compile_193 = Compile_106 + Compile_118 + \
      Compile_192 + Compile_30;
      std::complex<double> Compile_194 = -2*Compile_193*Compile_38;
      std::complex<double> Compile_197 = 6*lam1H3dUS*Compile_196;
      std::complex<double> Compile_198 = 7*Compile_176;
      std::complex<double> Compile_199 = Compile_192 + Compile_198;
      std::complex<double> Compile_200 = Compile_176*Compile_199;
      std::complex<double> Compile_201 = Compile_197 + Compile_200;
      std::complex<double> Compile_202 = Compile_201*Compile_5;
      std::complex<double> Compile_203 = Compile_194 + Compile_202;
      std::complex<double> Compile_204 = 2*Compile_203*Compile_6;
      std::complex<double> Compile_205 = pow(Compile_196,2);
      std::complex<double> Compile_206 = Compile_205*Compile_52;
      std::complex<double> Compile_207 = Compile_190 + Compile_191 + \
      Compile_204 + Compile_206 + Compile_28;
      std::complex<double> Compile_208 = sqrt(Compile_207);
      std::complex<double> Compile_211 = -(phi1*phi2*Compile_176);
      std::complex<double> Compile_212 = m12R3dUS + Compile_211;
      std::complex<double> Compile_217 = Compile_176*Compile_91*Compile_92;
      std::complex<double> Compile_218 = Compile_178 + Compile_181 + \
      Compile_208 + Compile_217 + Compile_79 + Compile_80;
      std::complex<double> Compile_177 = Compile_176*Compile_5;
      std::complex<double> Compile_179 = Compile_177 + Compile_178 + \
      Compile_80;
      std::complex<double> Compile_180 = 8*Compile_179;
      std::complex<double> Compile_213 = pow(Compile_212,-2);
      std::complex<double> Compile_214 = 3*lam1H3dUS*Compile_5;
      std::complex<double> Compile_215 = (Compile_176*Compile_6)/2.;
      std::complex<double> Compile_216 = m13dUS + Compile_214 + \
      Compile_215;
      std::complex<double> Compile_231 = 6*lam1H3dUS*Compile_5;
      std::complex<double> Compile_232 = -6*lam2H3dUS*Compile_6;
      std::complex<double> Compile_233 = \
      -(Compile_176*Compile_91*Compile_92);
      std::complex<double> Compile_234 = Compile_208 + Compile_21 + \
      Compile_22 + Compile_231 + Compile_232 + Compile_233;
      std::complex<double> Compile_222 = 1./Compile_212;
      std::complex<double> Compile_164 = -Compile_5;
      std::complex<double> Compile_165 = Compile_164 + Compile_6;
      

      return 2*(-0.020833333333333332*(Compile_2*pow(Compile_3*Compile_7,1.5)) - (Compile_2*pow((pow(g13dUS,2) + Compile_3)*Compile_7,1.5))/96.) \
      - (Compile_2*pow((Compile_110 + 8*(-m12R3dUS + \
      lam5H3dUS*phi1*phi2)*Compile_113*Compile_142 + \
      pow(Compile_142,2)*Compile_144*Compile_146)/(16 + \
      pow(abs(Compile_113*Compile_142),2)),1.5))/12. - \
      (Compile_2*pow((Compile_180 + \
      Compile_213*Compile_216*pow(Compile_218,2) - 8*(Compile_178 + \
      Compile_181 + Compile_208 + lam4H3dUS*Compile_5 + lam5H3dUS*Compile_5 \
      - lam4H3dUS*Compile_6 - lam5H3dUS*Compile_6 + Compile_79 + Compile_80 \
      + Compile_82 + Compile_84))/(16 + \
      pow(abs(Compile_218*Compile_222),2)),1.5))/12. - \
      (Compile_2*pow((Compile_180 + 8*Compile_234 + \
      Compile_213*Compile_216*pow(Compile_234,2))/(16 + \
      pow(abs(Compile_222*(Compile_165*Compile_176 + Compile_208 + \
      Compile_21 + Compile_22 + Compile_231 + Compile_232)),2)),1.5))/12. - \
      (Compile_2*pow((Compile_110 + 8*Compile_159 + \
      Compile_144*Compile_146*pow(Compile_159,2))/(16 + \
      pow(abs(Compile_113*(Compile_141 + Compile_107*Compile_165 + \
      Compile_21 + Compile_22 + Compile_23 + Compile_25)),2)),1.5))/12. - \
      (Compile_2*pow((Compile_18 + Compile_19 + Compile_20 + 2*Compile_56 + \
      ((Compile_21 + Compile_23 + \
      Compile_26)*pow(Compile_56,2)*Compile_61)/2.)/(4 + \
      pow(abs((Compile_21 + Compile_22 + Compile_25 + Compile_26 + \
      Compile_32 + Compile_55)*Compile_69),2)),1.5))/6. - \
      (Compile_2*pow((Compile_18 + Compile_19 + Compile_20 - 2*(Compile_55 \
      + Compile_79 + Compile_80 + Compile_81 + Compile_82 + Compile_83 + \
      Compile_84) + Compile_61*(m13dUS + (lam3H3dUS*Compile_6)/2. + \
      Compile_87)*pow(Compile_94,2))/(4 + \
      pow(abs(Compile_69*Compile_94),2)),1.5))/6.;
    }



    std::complex<double> V2( Eigen::VectorXd phi, double T, std::vector<std::complex<double>> par ) {

      std::complex<double> g13dUS = sqrt(par[0]);
      std::complex<double> g23dUS = sqrt(par[1]);
      std::complex<double> g33dUS = sqrt(par[2]);
      std::complex<double> lam1H3dUS = par[3];
      std::complex<double> lam2H3dUS = par[4];
      std::complex<double> lam3H3dUS = par[5];
      std::complex<double> lam4H3dUS = par[6];
      std::complex<double> lam5H3dUS = par[7];
      std::complex<double> lam6H3dUS = par[8];
      std::complex<double> lam7H3dUS = par[9];
      std::complex<double> m12R3dUS = par[10];
      std::complex<double> m13dUS = par[11];
      std::complex<double> m23dUS = par[12];
      std::complex<double> mu = par[13];
      std::complex<double> mu3US = par[14];
    
      std::complex<double> phi1(phi[0]/sqrt(T + 1e-15), 0);
      std::complex<double> phi2(phi[1]/sqrt(T + 1e-15), 0);

      std::complex<double> Compile_2 = pow(g13dUS,2);
      std::complex<double> Compile_5 = pow(g23dUS,2);
      std::complex<double> Compile_6 = Compile_2 + Compile_5;
      std::complex<double> Compile_7 = 1./Compile_6;
      std::complex<double> Compile_8 = pow(pi,-2);
      std::complex<double> Compile_9 = pow(phi1,2);
      std::complex<double> Compile_10 = pow(phi2,2);
      std::complex<double> Compile_11 = Compile_10 + Compile_9;
      std::complex<double> Compile_3 = pow(g23dUS,4);
      std::complex<double> Compile_28 = 2*m13dUS;
      std::complex<double> Compile_29 = -2*m23dUS;
      std::complex<double> Compile_21 = lam4H3dUS + lam5H3dUS;
      std::complex<double> Compile_37 = -lam3H3dUS;
      std::complex<double> Compile_36 = 2*lam1H3dUS;
      std::complex<double> Compile_38 = Compile_36 + Compile_37;
      std::complex<double> Compile_48 = -2*lam2H3dUS;
      std::complex<double> Compile_49 = lam3H3dUS + Compile_48;
      std::complex<double> Compile_20 = -2*m12R3dUS;
      std::complex<double> Compile_22 = phi1*phi2*Compile_21;
      std::complex<double> Compile_23 = Compile_20 + Compile_22;
      std::complex<double> Compile_30 = 2*lam1H3dUS*Compile_9;
      std::complex<double> Compile_33 = lam3H3dUS*Compile_10;
      std::complex<double> Compile_31 = -(lam3H3dUS*Compile_9);
      std::complex<double> Compile_32 = -2*lam2H3dUS*Compile_10;
      std::complex<double> Compile_34 = pow(m12R3dUS,2);
      std::complex<double> Compile_35 = 16*Compile_34;
      std::complex<double> Compile_39 = Compile_38*Compile_9;
      std::complex<double> Compile_40 = Compile_28 + Compile_29 + Compile_39;
      std::complex<double> Compile_41 = pow(Compile_40,2);
      std::complex<double> Compile_42 = -16*m12R3dUS*phi1*phi2*Compile_21;
      std::complex<double> Compile_43 = -m23dUS;
      std::complex<double> Compile_44 = m13dUS + Compile_43;
      std::complex<double> Compile_45 = 2*lam2H3dUS;
      std::complex<double> Compile_46 = Compile_37 + Compile_45;
      std::complex<double> Compile_47 = -2*Compile_44*Compile_46;
      std::complex<double> Compile_50 = Compile_38*Compile_49;
      std::complex<double> Compile_51 = pow(Compile_21,2);
      std::complex<double> Compile_52 = 2*Compile_51;
      std::complex<double> Compile_53 = Compile_50 + Compile_52;
      std::complex<double> Compile_54 = Compile_53*Compile_9;
      std::complex<double> Compile_55 = Compile_47 + Compile_54;
      std::complex<double> Compile_56 = 2*Compile_10*Compile_55;
      std::complex<double> Compile_57 = pow(Compile_49,2);
      std::complex<double> Compile_58 = pow(phi2,4);
      std::complex<double> Compile_59 = Compile_57*Compile_58;
      std::complex<double> Compile_60 = Compile_35 + Compile_41 + \
      Compile_42 + Compile_56 + Compile_59;
      std::complex<double> Compile_61 = sqrt(Compile_60);
      std::complex<double> Compile_62 = Compile_28 + Compile_29 + \
      Compile_30 + Compile_31 + Compile_32 + Compile_33 + Compile_61;
      std::complex<double> Compile_66 = pow(Compile_62,2);
      std::complex<double> Compile_24 = pow(Compile_23,-4);
      std::complex<double> Compile_25 = 4*m23dUS;
      std::complex<double> Compile_26 = 2*lam3H3dUS*Compile_9;
      std::complex<double> Compile_27 = 4*lam2H3dUS*Compile_10;
      std::complex<double> Compile_63 = 2*Compile_62;
      std::complex<double> Compile_64 = pow(Compile_23,-2);
      std::complex<double> Compile_65 = Compile_28 + Compile_30 + Compile_33;
      std::complex<double> Compile_67 = \
      (Compile_64*Compile_65*Compile_66)/2.;
      std::complex<double> Compile_68 = Compile_25 + Compile_26 + \
      Compile_27 + Compile_63 + Compile_67;
      std::complex<double> Compile_72 = pow(Compile_23,2);
      std::complex<double> Compile_89 = lam3H3dUS + lam4H3dUS;
      std::complex<double> Compile_77 = 2*m12R3dUS;
      std::complex<double> Compile_78 = -(phi1*phi2*Compile_21);
      std::complex<double> Compile_79 = Compile_77 + Compile_78;
      std::complex<double> Compile_80 = 1./Compile_79;
      std::complex<double> Compile_81 = Compile_28 + Compile_29 + \
      Compile_32 + Compile_33 + Compile_39 + Compile_61;
      std::complex<double> Compile_82 = Compile_80*Compile_81;
      std::complex<double> Compile_83 = abs(Compile_82);
      std::complex<double> Compile_84 = pow(Compile_83,2);
      std::complex<double> Compile_85 = 4 + Compile_84;
      std::complex<double> Compile_86 = pow(Compile_85,-3);
      std::complex<double> Compile_17 = Compile_11*Compile_6;
      std::complex<double> Compile_18 = sqrt(Compile_17);
      std::complex<double> Compile_114 = Compile_64*Compile_66;
      std::complex<double> Compile_115 = 4 + Compile_114;
      std::complex<double> Compile_116 = pow(Compile_115,2);
      std::complex<double> Compile_117 = 1./Compile_85;
      std::complex<double> Compile_118 = Compile_117*Compile_68;
      std::complex<double> Compile_119 = sqrt(Compile_118);
      std::complex<double> Compile_120 = pow(Compile_85,-2);
      std::complex<double> Compile_69 = pow(Compile_23,4);
      std::complex<double> Compile_70 = 48*lam2H3dUS*Compile_69;
      std::complex<double> Compile_71 = lam3H3dUS + lam4H3dUS + lam5H3dUS;
      std::complex<double> Compile_94 = -2*m13dUS;
      std::complex<double> Compile_95 = 2*m23dUS;
      std::complex<double> Compile_96 = -2*lam1H3dUS*Compile_9;
      std::complex<double> Compile_97 = lam3H3dUS*Compile_9;
      std::complex<double> Compile_98 = 2*lam2H3dUS*Compile_10;
      std::complex<double> Compile_99 = -(lam3H3dUS*Compile_10);
      std::complex<double> Compile_126 = Compile_61 + Compile_94 + \
      Compile_95 + Compile_96 + Compile_97 + Compile_98 + Compile_99;
      std::complex<double> Compile_136 = -phi2;
      std::complex<double> Compile_137 = phi1 + Compile_136;
      std::complex<double> Compile_138 = phi1 + phi2;
      std::complex<double> Compile_139 = lam3H3dUS*Compile_137*Compile_138;
      std::complex<double> Compile_140 = Compile_139 + Compile_61 + \
      Compile_94 + Compile_95 + Compile_96 + Compile_98;
      std::complex<double> Compile_132 = -2*Compile_126;
      std::complex<double> Compile_133 = lam1H3dUS*Compile_9;
      std::complex<double> Compile_134 = (lam3H3dUS*Compile_10)/2.;
      std::complex<double> Compile_135 = m13dUS + Compile_133 + \
      Compile_134;
      std::complex<double> Compile_141 = pow(Compile_140,2);
      std::complex<double> Compile_142 = \
      Compile_135*Compile_141*Compile_64;
      std::complex<double> Compile_143 = Compile_132 + Compile_142 + \
      Compile_25 + Compile_26 + Compile_27;
      std::complex<double> Compile_88 = 8*lam2H3dUS*Compile_72;
      std::complex<double> Compile_127 = pow(Compile_126,2);
      std::complex<double> Compile_100 = -Compile_61;
      std::complex<double> Compile_154 = Compile_100 + Compile_28 + \
      Compile_29 + Compile_30 + Compile_31 + Compile_32 + Compile_33;
      std::complex<double> Compile_103 = 2*lam5H3dUS*Compile_72;
      std::complex<double> Compile_144 = Compile_140*Compile_80;
      std::complex<double> Compile_145 = abs(Compile_144);
      std::complex<double> Compile_146 = pow(Compile_145,2);
      std::complex<double> Compile_147 = 4 + Compile_146;
      std::complex<double> Compile_148 = pow(Compile_147,-3);
      std::complex<double> Compile_111 = -Compile_5;
      std::complex<double> Compile_112 = Compile_111 + Compile_2;
      std::complex<double> Compile_113 = pow(Compile_112,2);
      std::complex<double> Compile_122 = -Compile_2;
      std::complex<double> Compile_123 = Compile_122 + Compile_5;
      std::complex<double> Compile_124 = pow(Compile_123,2);
      std::complex<double> Compile_163 = Compile_141*Compile_64;
      std::complex<double> Compile_164 = 4 + Compile_163;
      std::complex<double> Compile_165 = pow(Compile_164,2);
      std::complex<double> Compile_166 = 1./Compile_147;
      std::complex<double> Compile_167 = Compile_143*Compile_166;
      std::complex<double> Compile_168 = sqrt(Compile_167);
      std::complex<double> Compile_169 = pow(Compile_147,-2);
      std::complex<double> Compile_172 = pow(m23dUS,2);
      std::complex<double> Compile_184 = pow(lam3H3dUS,2);
      std::complex<double> Compile_191 = pow(phi1,4);
      std::complex<double> Compile_190 = pow(lam1H3dUS,2);
      std::complex<double> Compile_175 = -3*lam5H3dUS;
      std::complex<double> Compile_176 = lam1H3dUS + lam2H3dUS + lam3H3dUS \
      + lam4H3dUS + Compile_175;
      std::complex<double> Compile_194 = pow(lam3H3dUS,3);
      std::complex<double> Compile_217 = pow(lam4H3dUS,2);
      std::complex<double> Compile_231 = pow(lam5H3dUS,2);
      std::complex<double> Compile_238 = pow(lam2H3dUS,2);
      std::complex<double> Compile_179 = -lam5H3dUS;
      std::complex<double> Compile_180 = lam3H3dUS + lam4H3dUS + Compile_179;
      std::complex<double> Compile_173 = 4*lam3H3dUS*Compile_172;
      std::complex<double> Compile_174 = 4*lam4H3dUS*Compile_172;
      std::complex<double> Compile_178 = pow(m13dUS,2);
      std::complex<double> Compile_183 = \
      -8*m23dUS*lam1H3dUS*lam3H3dUS*Compile_9;
      std::complex<double> Compile_185 = 4*m23dUS*Compile_184*Compile_9;
      std::complex<double> Compile_186 = \
      -8*m23dUS*lam1H3dUS*lam4H3dUS*Compile_9;
      std::complex<double> Compile_187 = \
      4*m23dUS*lam3H3dUS*lam4H3dUS*Compile_9;
      std::complex<double> Compile_192 = 4*lam3H3dUS*Compile_190*Compile_191;
      std::complex<double> Compile_193 = \
      -4*lam1H3dUS*Compile_184*Compile_191;
      std::complex<double> Compile_195 = Compile_191*Compile_194;
      std::complex<double> Compile_196 = 4*lam4H3dUS*Compile_190*Compile_191;
      std::complex<double> Compile_197 = \
      -4*lam1H3dUS*lam3H3dUS*lam4H3dUS*Compile_191;
      std::complex<double> Compile_198 = lam4H3dUS*Compile_184*Compile_191;
      std::complex<double> Compile_252 = 3*lam1H3dUS;
      std::complex<double> Compile_253 = 3*lam2H3dUS;
      std::complex<double> Compile_254 = -lam4H3dUS;
      std::complex<double> Compile_255 = Compile_179 + Compile_252 + \
      Compile_253 + Compile_254 + Compile_37;
      std::complex<double> Compile_203 = \
      8*m23dUS*lam2H3dUS*lam3H3dUS*Compile_10;
      std::complex<double> Compile_204 = -4*m23dUS*Compile_10*Compile_184;
      std::complex<double> Compile_205 = \
      8*m23dUS*lam2H3dUS*lam4H3dUS*Compile_10;
      std::complex<double> Compile_206 = \
      -4*m23dUS*lam3H3dUS*lam4H3dUS*Compile_10;
      std::complex<double> Compile_209 = \
      -8*lam1H3dUS*lam2H3dUS*lam3H3dUS*Compile_10*Compile_9;
      std::complex<double> Compile_210 = \
      4*lam1H3dUS*Compile_10*Compile_184*Compile_9;
      std::complex<double> Compile_211 = \
      4*lam2H3dUS*Compile_10*Compile_184*Compile_9;
      std::complex<double> Compile_212 = -2*Compile_10*Compile_194*Compile_9;
      std::complex<double> Compile_213 = \
      -8*lam1H3dUS*lam2H3dUS*lam4H3dUS*Compile_10*Compile_9;
      std::complex<double> Compile_214 = \
      4*lam1H3dUS*lam3H3dUS*lam4H3dUS*Compile_10*Compile_9;
      std::complex<double> Compile_215 = \
      4*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_10*Compile_9;
      std::complex<double> Compile_216 = \
      -2*lam4H3dUS*Compile_10*Compile_184*Compile_9;
      std::complex<double> Compile_221 = pow(lam4H3dUS,3);
      std::complex<double> Compile_236 = pow(lam5H3dUS,3);
      std::complex<double> Compile_239 = \
      4*lam3H3dUS*Compile_238*Compile_58;
      std::complex<double> Compile_240 = -4*lam2H3dUS*Compile_184*Compile_58;
      std::complex<double> Compile_241 = Compile_194*Compile_58;
      std::complex<double> Compile_242 = \
      4*lam4H3dUS*Compile_238*Compile_58;
      std::complex<double> Compile_243 = \
      -4*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_58;
      std::complex<double> Compile_244 = lam4H3dUS*Compile_184*Compile_58;
      std::complex<double> Compile_248 = Compile_95 + Compile_96 + \
      Compile_97 + Compile_98 + Compile_99;
      std::complex<double> Compile_290 = -(lam5H3dUS*phi1*phi2);
      std::complex<double> Compile_291 = m12R3dUS + Compile_290;
      std::complex<double> Compile_300 = Compile_137*Compile_138*Compile_180;
      std::complex<double> Compile_301 = lam5H3dUS + Compile_254 + \
      Compile_36 + Compile_37;
      std::complex<double> Compile_302 = Compile_301*Compile_9;
      std::complex<double> Compile_303 = Compile_28 + Compile_29 + \
      Compile_302;
      std::complex<double> Compile_304 = pow(Compile_303,2);
      std::complex<double> Compile_305 = -32*m12R3dUS*lam5H3dUS*phi1*phi2;
      std::complex<double> Compile_306 = lam5H3dUS + Compile_254 + \
      Compile_37 + Compile_45;
      std::complex<double> Compile_307 = 2*Compile_306*Compile_44;
      std::complex<double> Compile_308 = Compile_254 + Compile_36 + \
      Compile_37;
      std::complex<double> Compile_309 = lam3H3dUS + lam4H3dUS + \
      Compile_48;
      std::complex<double> Compile_310 = -(Compile_308*Compile_309);
      std::complex<double> Compile_311 = lam1H3dUS + lam2H3dUS + \
      Compile_254 + Compile_37;
      std::complex<double> Compile_312 = 2*lam5H3dUS*Compile_311;
      std::complex<double> Compile_313 = -7*Compile_231;
      std::complex<double> Compile_314 = Compile_310 + Compile_312 + \
      Compile_313;
      std::complex<double> Compile_315 = Compile_314*Compile_9;
      std::complex<double> Compile_316 = Compile_307 + Compile_315;
      std::complex<double> Compile_317 = -2*Compile_10*Compile_316;
      std::complex<double> Compile_318 = lam3H3dUS + lam4H3dUS + \
      Compile_179 + Compile_48;
      std::complex<double> Compile_319 = pow(Compile_318,2);
      std::complex<double> Compile_320 = Compile_319*Compile_58;
      std::complex<double> Compile_321 = Compile_304 + Compile_305 + \
      Compile_317 + Compile_320 + Compile_35;
      std::complex<double> Compile_322 = sqrt(Compile_321);
      std::complex<double> Compile_323 = Compile_300 + Compile_322 + \
      Compile_94 + Compile_95 + Compile_96 + Compile_98;
      std::complex<double> Compile_328 = pow(Compile_323,2);
      std::complex<double> Compile_296 = 1./Compile_291;
      std::complex<double> Compile_15 = Compile_11*Compile_5;
      std::complex<double> Compile_16 = sqrt(Compile_15);
      std::complex<double> Compile_338 = Compile_296*Compile_323;
      std::complex<double> Compile_339 = abs(Compile_338);
      std::complex<double> Compile_340 = pow(Compile_339,2);
      std::complex<double> Compile_341 = 16 + Compile_340;
      std::complex<double> Compile_344 = Compile_100 + Compile_139 + \
      Compile_94 + Compile_95 + Compile_96 + Compile_98;
      std::complex<double> Compile_348 = 1./Compile_341;
      std::complex<double> Compile_345 = \
      -(Compile_296*Compile_323*Compile_344*Compile_80);
      std::complex<double> Compile_346 = -8 + Compile_345;
      std::complex<double> Compile_347 = pow(Compile_346,2);
      std::complex<double> Compile_293 = Compile_180*Compile_9;
      std::complex<double> Compile_294 = Compile_293 + Compile_95 + \
      Compile_98;
      std::complex<double> Compile_295 = 8*Compile_294;
      std::complex<double> Compile_297 = -m12R3dUS;
      std::complex<double> Compile_298 = lam5H3dUS*phi1*phi2;
      std::complex<double> Compile_299 = Compile_297 + Compile_298;
      std::complex<double> Compile_324 = \
      8*Compile_296*Compile_299*Compile_323;
      std::complex<double> Compile_325 = pow(Compile_291,-2);
      std::complex<double> Compile_326 = (Compile_10*Compile_180)/2.;
      std::complex<double> Compile_327 = m13dUS + Compile_133 + \
      Compile_326;
      std::complex<double> Compile_329 = Compile_325*Compile_327*Compile_328;
      std::complex<double> Compile_330 = Compile_295 + Compile_324 + \
      Compile_329;
      std::complex<double> Compile_350 = \
      Compile_296*Compile_323*Compile_344*Compile_80;
      std::complex<double> Compile_351 = 8 + Compile_350;
      std::complex<double> Compile_352 = pow(Compile_351,2);
      std::complex<double> Compile_362 = Compile_330*Compile_348;
      std::complex<double> Compile_363 = sqrt(Compile_362);
      std::complex<double> Compile_367 = \
      -(Compile_137*Compile_138*Compile_180);
      std::complex<double> Compile_368 = -Compile_322;
      std::complex<double> Compile_369 = Compile_28 + Compile_29 + \
      Compile_30 + Compile_32 + Compile_367 + Compile_368;
      std::complex<double> Compile_379 = Compile_21*Compile_291*Compile_62;
      std::complex<double> Compile_380 = lam3H3dUS*Compile_369*Compile_79;
      std::complex<double> Compile_381 = Compile_379 + Compile_380;
      std::complex<double> Compile_373 = 2*Compile_21*Compile_291*Compile_79;
      std::complex<double> Compile_374 = lam1H3dUS*Compile_369*Compile_62;
      std::complex<double> Compile_375 = Compile_373 + Compile_374;
      std::complex<double> Compile_366 = \
      4*lam3H3dUS*Compile_291*Compile_62;
      std::complex<double> Compile_370 = Compile_21*Compile_369*Compile_79;
      std::complex<double> Compile_371 = Compile_366 + Compile_370;
      std::complex<double> Compile_383 = 32*lam2H3dUS*Compile_291*Compile_79;
      std::complex<double> Compile_384 = Compile_21*Compile_369*Compile_62;
      std::complex<double> Compile_385 = Compile_383 + Compile_384;
      std::complex<double> Compile_354 = \
      -(Compile_140*Compile_296*Compile_323*Compile_80);
      std::complex<double> Compile_355 = -8 + Compile_354;
      std::complex<double> Compile_356 = pow(Compile_355,2);
      std::complex<double> Compile_358 = \
      Compile_140*Compile_296*Compile_323*Compile_80;
      std::complex<double> Compile_359 = 8 + Compile_358;
      std::complex<double> Compile_360 = pow(Compile_359,2);
      std::complex<double> Compile_411 = \
      Compile_154*Compile_21*Compile_291;
      std::complex<double> Compile_412 = Compile_380 + Compile_411;
      std::complex<double> Compile_406 = lam1H3dUS*Compile_154*Compile_369;
      std::complex<double> Compile_407 = Compile_373 + Compile_406;
      std::complex<double> Compile_403 = \
      -4*lam3H3dUS*Compile_126*Compile_291;
      std::complex<double> Compile_404 = Compile_370 + Compile_403;
      std::complex<double> Compile_414 = \
      Compile_154*Compile_21*Compile_369;
      std::complex<double> Compile_415 = Compile_383 + Compile_414;
      std::complex<double> Compile_292 = pow(Compile_291,-4);
      std::complex<double> Compile_431 = Compile_28 + Compile_29 + \
      Compile_30 + Compile_32 + Compile_322 + Compile_367;
      std::complex<double> Compile_331 = pow(Compile_291,4);
      std::complex<double> Compile_332 = 768*lam2H3dUS*Compile_331;
      std::complex<double> Compile_333 = pow(Compile_291,2);
      std::complex<double> Compile_433 = pow(Compile_431,2);
      std::complex<double> Compile_440 = -Compile_9;
      std::complex<double> Compile_441 = Compile_10 + Compile_440;
      std::complex<double> Compile_442 = Compile_180*Compile_441;
      std::complex<double> Compile_443 = Compile_28 + Compile_29 + \
      Compile_30 + Compile_32 + Compile_322 + Compile_442;
      std::complex<double> Compile_444 = Compile_296*Compile_443;
      std::complex<double> Compile_445 = abs(Compile_444);
      std::complex<double> Compile_446 = pow(Compile_445,2);
      std::complex<double> Compile_447 = 16 + Compile_446;
      std::complex<double> Compile_450 = Compile_300 + Compile_368 + \
      Compile_94 + Compile_95 + Compile_96 + Compile_98;
      std::complex<double> Compile_454 = 1./Compile_447;
      std::complex<double> Compile_451 = \
      -(Compile_296*Compile_344*Compile_450*Compile_80);
      std::complex<double> Compile_452 = -8 + Compile_451;
      std::complex<double> Compile_453 = pow(Compile_452,2);
      std::complex<double> Compile_432 = 8*Compile_431;
      std::complex<double> Compile_434 = Compile_325*Compile_327*Compile_433;
      std::complex<double> Compile_435 = Compile_295 + Compile_432 + \
      Compile_434;
      std::complex<double> Compile_456 = \
      Compile_296*Compile_344*Compile_450*Compile_80;
      std::complex<double> Compile_457 = 8 + Compile_456;
      std::complex<double> Compile_458 = pow(Compile_457,2);
      std::complex<double> Compile_468 = Compile_435*Compile_454;
      std::complex<double> Compile_469 = sqrt(Compile_468);
      std::complex<double> Compile_480 = lam3H3dUS*Compile_431*Compile_79;
      std::complex<double> Compile_481 = Compile_379 + Compile_480;
      std::complex<double> Compile_475 = lam1H3dUS*Compile_431*Compile_62;
      std::complex<double> Compile_476 = Compile_373 + Compile_475;
      std::complex<double> Compile_472 = Compile_21*Compile_431*Compile_79;
      std::complex<double> Compile_473 = Compile_366 + Compile_472;
      std::complex<double> Compile_483 = Compile_21*Compile_431*Compile_62;
      std::complex<double> Compile_484 = Compile_383 + Compile_483;
      std::complex<double> Compile_460 = \
      -(Compile_140*Compile_296*Compile_450*Compile_80);
      std::complex<double> Compile_461 = -8 + Compile_460;
      std::complex<double> Compile_462 = pow(Compile_461,2);
      std::complex<double> Compile_464 = \
      Compile_140*Compile_296*Compile_450*Compile_80;
      std::complex<double> Compile_465 = 8 + Compile_464;
      std::complex<double> Compile_466 = pow(Compile_465,2);
      std::complex<double> Compile_509 = Compile_411 + Compile_480;
      std::complex<double> Compile_504 = lam1H3dUS*Compile_154*Compile_431;
      std::complex<double> Compile_505 = Compile_373 + Compile_504;
      std::complex<double> Compile_502 = Compile_403 + Compile_472;
      std::complex<double> Compile_511 = \
      Compile_154*Compile_21*Compile_431;
      std::complex<double> Compile_512 = Compile_383 + Compile_511;
      std::complex<double> Compile_256 = 8*Compile_255*Compile_34;
      std::complex<double> Compile_257 = 4*lam5H3dUS*Compile_172;
      std::complex<double> Compile_258 = 4*Compile_178*Compile_71;
      std::complex<double> Compile_259 = \
      -8*m23dUS*lam1H3dUS*lam5H3dUS*Compile_9;
      std::complex<double> Compile_261 = 4*lam5H3dUS*Compile_190*Compile_191;
      std::complex<double> Compile_201 = \
      -(lam5H3dUS*Compile_184*Compile_191);
      std::complex<double> Compile_265 = \
      8*m23dUS*lam2H3dUS*lam5H3dUS*Compile_10;
      std::complex<double> Compile_270 = -2*Compile_10*Compile_221*Compile_9;
      std::complex<double> Compile_271 = \
      -8*lam1H3dUS*lam2H3dUS*lam5H3dUS*Compile_10*Compile_9;
      std::complex<double> Compile_226 = \
      2*lam5H3dUS*Compile_10*Compile_184*Compile_9;
      std::complex<double> Compile_229 = \
      4*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_10*Compile_9;
      std::complex<double> Compile_282 = \
      -6*lam4H3dUS*Compile_10*Compile_231*Compile_9;
      std::complex<double> Compile_284 = \
      4*lam5H3dUS*Compile_238*Compile_58;
      std::complex<double> Compile_247 = -(lam5H3dUS*Compile_184*Compile_58);
      std::complex<double> Compile_568 = lam4H3dUS*Compile_9;
      std::complex<double> Compile_580 = 6*lam2H3dUS*Compile_10;
      std::complex<double> Compile_570 = -(lam4H3dUS*Compile_10);
      std::complex<double> Compile_592 = 6*lam2H3dUS;
      std::complex<double> Compile_595 = -6*lam2H3dUS;
      std::complex<double> Compile_596 = lam3H3dUS + lam4H3dUS + lam5H3dUS \
      + Compile_595;
      std::complex<double> Compile_576 = -(phi1*phi2*Compile_71);
      std::complex<double> Compile_577 = m12R3dUS + Compile_576;
      std::complex<double> Compile_583 = -6*lam1H3dUS*Compile_9;
      std::complex<double> Compile_586 = 6*lam1H3dUS;
      std::complex<double> Compile_587 = Compile_179 + Compile_254 + \
      Compile_37 + Compile_586;
      std::complex<double> Compile_588 = Compile_587*Compile_9;
      std::complex<double> Compile_589 = Compile_28 + Compile_29 + \
      Compile_588;
      std::complex<double> Compile_590 = pow(Compile_589,2);
      std::complex<double> Compile_591 = -32*m12R3dUS*phi1*phi2*Compile_71;
      std::complex<double> Compile_593 = Compile_179 + Compile_254 + \
      Compile_37 + Compile_592;
      std::complex<double> Compile_594 = -2*Compile_44*Compile_593;
      std::complex<double> Compile_597 = 6*lam1H3dUS*Compile_596;
      std::complex<double> Compile_598 = 7*Compile_71;
      std::complex<double> Compile_599 = Compile_592 + Compile_598;
      std::complex<double> Compile_600 = Compile_599*Compile_71;
      std::complex<double> Compile_601 = Compile_597 + Compile_600;
      std::complex<double> Compile_602 = Compile_601*Compile_9;
      std::complex<double> Compile_603 = Compile_594 + Compile_602;
      std::complex<double> Compile_604 = 2*Compile_10*Compile_603;
      std::complex<double> Compile_605 = pow(Compile_596,2);
      std::complex<double> Compile_606 = Compile_58*Compile_605;
      std::complex<double> Compile_607 = Compile_35 + Compile_590 + \
      Compile_591 + Compile_604 + Compile_606;
      std::complex<double> Compile_608 = sqrt(Compile_607);
      std::complex<double> Compile_615 = \
      Compile_137*Compile_138*Compile_71;
      std::complex<double> Compile_616 = Compile_580 + Compile_583 + \
      Compile_608 + Compile_615 + Compile_94 + Compile_95;
      std::complex<double> Compile_617 = pow(Compile_616,2);
      std::complex<double> Compile_627 = 1./Compile_577;
      std::complex<double> Compile_628 = Compile_616*Compile_627;
      std::complex<double> Compile_629 = abs(Compile_628);
      std::complex<double> Compile_630 = pow(Compile_629,2);
      std::complex<double> Compile_631 = 16 + Compile_630;
      std::complex<double> Compile_637 = 1./Compile_631;
      std::complex<double> Compile_579 = Compile_71*Compile_9;
      std::complex<double> Compile_581 = Compile_579 + Compile_580 + \
      Compile_95;
      std::complex<double> Compile_582 = 8*Compile_581;
      std::complex<double> Compile_584 = lam5H3dUS*Compile_9;
      std::complex<double> Compile_585 = -(lam5H3dUS*Compile_10);
      std::complex<double> Compile_609 = Compile_568 + Compile_570 + \
      Compile_580 + Compile_583 + Compile_584 + Compile_585 + Compile_608 + \
      Compile_94 + Compile_95 + Compile_97 + Compile_99;
      std::complex<double> Compile_610 = -8*Compile_609;
      std::complex<double> Compile_611 = pow(Compile_577,-2);
      std::complex<double> Compile_612 = 3*lam1H3dUS*Compile_9;
      std::complex<double> Compile_613 = (Compile_10*Compile_71)/2.;
      std::complex<double> Compile_614 = m13dUS + Compile_612 + \
      Compile_613;
      std::complex<double> Compile_618 = Compile_611*Compile_614*Compile_617;
      std::complex<double> Compile_619 = Compile_582 + Compile_610 + \
      Compile_618;
      std::complex<double> Compile_657 = 6*lam1H3dUS*Compile_9;
      std::complex<double> Compile_658 = -6*lam2H3dUS*Compile_10;
      std::complex<double> Compile_659 = \
      -(Compile_137*Compile_138*Compile_71);
      std::complex<double> Compile_660 = -Compile_608;
      std::complex<double> Compile_661 = Compile_28 + Compile_29 + \
      Compile_657 + Compile_658 + Compile_659 + Compile_660;
      std::complex<double> Compile_654 = Compile_619*Compile_637;
      std::complex<double> Compile_655 = sqrt(Compile_654);
      std::complex<double> Compile_672 = \
      4*lam3H3dUS*Compile_577*Compile_62;
      std::complex<double> Compile_673 = Compile_21*Compile_661*Compile_79;
      std::complex<double> Compile_674 = Compile_672 + Compile_673;
      std::complex<double> Compile_666 = 2*Compile_21*Compile_577*Compile_79;
      std::complex<double> Compile_667 = lam1H3dUS*Compile_62*Compile_661;
      std::complex<double> Compile_668 = Compile_666 + Compile_667;
      std::complex<double> Compile_662 = Compile_21*Compile_577*Compile_62;
      std::complex<double> Compile_663 = lam3H3dUS*Compile_661*Compile_79;
      std::complex<double> Compile_664 = Compile_662 + Compile_663;
      std::complex<double> Compile_676 = 32*lam2H3dUS*Compile_577*Compile_79;
      std::complex<double> Compile_677 = Compile_21*Compile_62*Compile_661;
      std::complex<double> Compile_678 = Compile_676 + Compile_677;
      std::complex<double> Compile_706 = 4*lam3H3dUS*Compile_154*Compile_577;
      std::complex<double> Compile_707 = Compile_673 + Compile_706;
      std::complex<double> Compile_701 = lam1H3dUS*Compile_154*Compile_661;
      std::complex<double> Compile_702 = Compile_666 + Compile_701;
      std::complex<double> Compile_698 = \
      Compile_154*Compile_21*Compile_577;
      std::complex<double> Compile_699 = Compile_663 + Compile_698;
      std::complex<double> Compile_709 = \
      Compile_154*Compile_21*Compile_661;
      std::complex<double> Compile_710 = Compile_676 + Compile_709;
      std::complex<double> Compile_643 = \
      Compile_296*Compile_323*Compile_616*Compile_627;
      std::complex<double> Compile_644 = 16 + Compile_643;
      std::complex<double> Compile_645 = pow(Compile_644,2);
      std::complex<double> Compile_737 = Compile_180*Compile_369*Compile_577;
      std::complex<double> Compile_738 = lam5H3dUS*Compile_291*Compile_661;
      std::complex<double> Compile_739 = Compile_737 + Compile_738;
      std::complex<double> Compile_731 = 8*lam5H3dUS*Compile_291*Compile_577;
      std::complex<double> Compile_732 = lam1H3dUS*Compile_369*Compile_661;
      std::complex<double> Compile_733 = Compile_731 + Compile_732;
      std::complex<double> Compile_727 = lam5H3dUS*Compile_369*Compile_577;
      std::complex<double> Compile_728 = Compile_180*Compile_291*Compile_661;
      std::complex<double> Compile_729 = Compile_727 + Compile_728;
      std::complex<double> Compile_741 = \
      32*lam2H3dUS*Compile_291*Compile_577;
      std::complex<double> Compile_742 = lam5H3dUS*Compile_369*Compile_661;
      std::complex<double> Compile_743 = Compile_741 + Compile_742;
      std::complex<double> Compile_647 = \
      Compile_296*Compile_450*Compile_616*Compile_627;
      std::complex<double> Compile_648 = 16 + Compile_647;
      std::complex<double> Compile_649 = pow(Compile_648,2);
      std::complex<double> Compile_768 = Compile_180*Compile_431*Compile_577;
      std::complex<double> Compile_769 = Compile_738 + Compile_768;
      std::complex<double> Compile_763 = lam1H3dUS*Compile_431*Compile_661;
      std::complex<double> Compile_764 = Compile_731 + Compile_763;
      std::complex<double> Compile_760 = lam5H3dUS*Compile_431*Compile_577;
      std::complex<double> Compile_761 = Compile_728 + Compile_760;
      std::complex<double> Compile_771 = lam5H3dUS*Compile_431*Compile_661;
      std::complex<double> Compile_772 = Compile_741 + Compile_771;
      std::complex<double> Compile_578 = pow(Compile_577,-4);
      std::complex<double> Compile_788 = Compile_28 + Compile_29 + \
      Compile_608 + Compile_657 + Compile_658 + Compile_659;
      std::complex<double> Compile_620 = pow(Compile_577,4);
      std::complex<double> Compile_621 = 768*lam2H3dUS*Compile_620;
      std::complex<double> Compile_622 = pow(Compile_577,2);
      std::complex<double> Compile_790 = pow(Compile_788,2);
      std::complex<double> Compile_797 = Compile_441*Compile_71;
      std::complex<double> Compile_798 = Compile_28 + Compile_29 + \
      Compile_608 + Compile_657 + Compile_658 + Compile_797;
      std::complex<double> Compile_799 = Compile_627*Compile_798;
      std::complex<double> Compile_800 = abs(Compile_799);
      std::complex<double> Compile_801 = pow(Compile_800,2);
      std::complex<double> Compile_802 = 16 + Compile_801;
      std::complex<double> Compile_805 = Compile_580 + Compile_583 + \
      Compile_615 + Compile_660 + Compile_94 + Compile_95;
      std::complex<double> Compile_809 = 1./Compile_802;
      std::complex<double> Compile_789 = 8*Compile_788;
      std::complex<double> Compile_791 = Compile_611*Compile_614*Compile_790;
      std::complex<double> Compile_792 = Compile_582 + Compile_789 + \
      Compile_791;
      std::complex<double> Compile_826 = Compile_792*Compile_809;
      std::complex<double> Compile_827 = sqrt(Compile_826);
      std::complex<double> Compile_837 = Compile_21*Compile_788*Compile_79;
      std::complex<double> Compile_838 = Compile_672 + Compile_837;
      std::complex<double> Compile_832 = lam1H3dUS*Compile_62*Compile_788;
      std::complex<double> Compile_833 = Compile_666 + Compile_832;
      std::complex<double> Compile_829 = lam3H3dUS*Compile_788*Compile_79;
      std::complex<double> Compile_830 = Compile_662 + Compile_829;
      std::complex<double> Compile_840 = Compile_21*Compile_62*Compile_788;
      std::complex<double> Compile_841 = Compile_676 + Compile_840;
      std::complex<double> Compile_868 = Compile_706 + Compile_837;
      std::complex<double> Compile_863 = lam1H3dUS*Compile_154*Compile_788;
      std::complex<double> Compile_864 = Compile_666 + Compile_863;
      std::complex<double> Compile_861 = Compile_698 + Compile_829;
      std::complex<double> Compile_870 = \
      Compile_154*Compile_21*Compile_788;
      std::complex<double> Compile_871 = Compile_676 + Compile_870;
      std::complex<double> Compile_815 = \
      Compile_296*Compile_323*Compile_627*Compile_805;
      std::complex<double> Compile_816 = 16 + Compile_815;
      std::complex<double> Compile_817 = pow(Compile_816,2);
      std::complex<double> Compile_896 = lam5H3dUS*Compile_291*Compile_788;
      std::complex<double> Compile_897 = Compile_737 + Compile_896;
      std::complex<double> Compile_891 = lam1H3dUS*Compile_369*Compile_788;
      std::complex<double> Compile_892 = Compile_731 + Compile_891;
      std::complex<double> Compile_888 = Compile_180*Compile_291*Compile_788;
      std::complex<double> Compile_889 = Compile_727 + Compile_888;
      std::complex<double> Compile_899 = lam5H3dUS*Compile_369*Compile_788;
      std::complex<double> Compile_900 = Compile_741 + Compile_899;
      std::complex<double> Compile_819 = \
      Compile_296*Compile_450*Compile_627*Compile_805;
      std::complex<double> Compile_820 = 16 + Compile_819;
      std::complex<double> Compile_821 = pow(Compile_820,2);
      std::complex<double> Compile_924 = Compile_768 + Compile_896;
      std::complex<double> Compile_919 = lam1H3dUS*Compile_431*Compile_788;
      std::complex<double> Compile_920 = Compile_731 + Compile_919;
      std::complex<double> Compile_917 = Compile_760 + Compile_888;
      std::complex<double> Compile_926 = lam5H3dUS*Compile_431*Compile_788;
      std::complex<double> Compile_927 = Compile_741 + Compile_926;
      std::complex<double> Compile_569 = -(lam5H3dUS*Compile_9);
      std::complex<double> Compile_571 = lam5H3dUS*Compile_10;
      std::complex<double> Compile_13 = pow(g23dUS,6);
      std::complex<double> Compile_989 = 1./sqrt(Compile_15);
      std::complex<double> Compile_990 = 2*mu3US*Compile_989;
      std::complex<double> Compile_991 = log(Compile_990);
      std::complex<double> Compile_992 = 0.5 + Compile_991;
      std::complex<double> Compile_997 = pow(Compile_11,3);
      std::complex<double> Compile_1013 = pow(Compile_11,4);
      std::complex<double> Compile_1019 = pow(Compile_11,2);
      std::complex<double> Compile_1015 = pow(Compile_6,2);
      std::complex<double> Compile_1031 = Compile_8/32.;
      std::complex<double> Compile_1006 = 1./sqrt(Compile_17);
      std::complex<double> Compile_1007 = 2*mu3US*Compile_1006;
      std::complex<double> Compile_1008 = log(Compile_1007);
      std::complex<double> Compile_1009 = 0.5 + Compile_1008;
      std::complex<double> Compile_1044 = (Compile_1019*Compile_3)/16.;
      std::complex<double> Compile_1045 = \
      (3*Compile_1019*Compile_5*Compile_6)/8.;
      std::complex<double> Compile_1046 = (Compile_1015*Compile_1019)/16.;
      std::complex<double> Compile_1047 = Compile_1044 + Compile_1045 + \
      Compile_1046;
      std::complex<double> Compile_1067 = (Compile_1019*Compile_3)/2.;
      std::complex<double> Compile_1068 = \
      (3*Compile_1019*Compile_5*Compile_6)/4.;
      std::complex<double> Compile_1069 = Compile_1046 + Compile_1067 + \
      Compile_1068;
      std::complex<double> Compile_1052 = Compile_18/2.;
      std::complex<double> Compile_1011 = pow(Compile_6,-2);
      std::complex<double> Compile_1012 = pow(Compile_11,-3);
      std::complex<double> Compile_1014 = \
      (Compile_1013*Compile_13*Compile_6*Compile_8)/8192.;
      std::complex<double> Compile_1016 = \
      (Compile_1013*Compile_1015*Compile_3*Compile_8)/16384.;
      std::complex<double> Compile_1017 = Compile_1014 + Compile_1016;
      std::complex<double> Compile_1018 = (40*Compile_1017)/3.;
      std::complex<double> Compile_1020 = (3*Compile_1019*Compile_3)/16.;
      std::complex<double> Compile_1021 = \
      (-11*Compile_1019*Compile_5*Compile_6)/16.;
      std::complex<double> Compile_1022 = \
      (-15*Compile_1015*Compile_1019)/16.;
      std::complex<double> Compile_1023 = Compile_1020 + Compile_1021 + \
      Compile_1022;
      std::complex<double> Compile_1024 = \
      (Compile_1023*Compile_11*Compile_16*Compile_18*Compile_5*Compile_8)/\
      128.;
      std::complex<double> Compile_1025 = (-7*Compile_1019*Compile_3)/2.;
      std::complex<double> Compile_1026 = \
      (15*Compile_1019*Compile_5*Compile_6)/8.;
      std::complex<double> Compile_1027 = \
      (3*Compile_1015*Compile_1019)/16.;
      std::complex<double> Compile_1028 = Compile_1025 + Compile_1026 + \
      Compile_1027;
      std::complex<double> Compile_1029 = \
      (Compile_1019*Compile_1028*Compile_5*Compile_6*Compile_8)/256.;
      std::complex<double> Compile_1030 = pow(g23dUS,8);
      std::complex<double> Compile_1032 = (-3*Compile_8*Compile_992)/16.;
      std::complex<double> Compile_1033 = Compile_1031 + Compile_1032;
      std::complex<double> Compile_1034 = \
      (Compile_1013*Compile_1030*Compile_1033)/128.;
      std::complex<double> Compile_1035 = pow(Compile_6,4);
      std::complex<double> Compile_1036 = (-3*Compile_1009*Compile_8)/16.;
      std::complex<double> Compile_1037 = Compile_1031 + Compile_1036;
      std::complex<double> Compile_1038 = \
      (Compile_1013*Compile_1035*Compile_1037)/256.;
      std::complex<double> Compile_1043 = \
      (-3*Compile_1019*Compile_5*Compile_6)/2.;
      std::complex<double> Compile_1048 = -2*Compile_1047;
      std::complex<double> Compile_1049 = Compile_1043 + Compile_1048;
      std::complex<double> Compile_1050 = (Compile_1049*Compile_8)/64.;
      std::complex<double> Compile_1051 = Compile_16/2.;
      std::complex<double> Compile_1053 = Compile_1051 + Compile_1052;
      std::complex<double> Compile_1054 = 1./Compile_1053;
      std::complex<double> Compile_1055 = mu3US*Compile_1054;
      std::complex<double> Compile_1056 = log(Compile_1055);
      std::complex<double> Compile_1057 = 0.5 + Compile_1056;
      std::complex<double> Compile_1058 = \
      (3*Compile_1047*Compile_1057*Compile_8)/16.;
      std::complex<double> Compile_1059 = Compile_1050 + Compile_1058;
      std::complex<double> Compile_1085 = (Compile_11*Compile_5)/4.;
      std::complex<double> Compile_1086 = -0.25*(Compile_11*Compile_6);
      std::complex<double> Compile_1087 = Compile_1085 + Compile_1086;
      std::complex<double> Compile_1088 = pow(Compile_1087,2);
      std::complex<double> Compile_1040 = (Compile_11*Compile_6)/4.;
      std::complex<double> Compile_1090 = Compile_1040 + Compile_1085;
      std::complex<double> Compile_1097 = (Compile_1019*Compile_3)/8.;
      std::complex<double> Compile_1098 = \
      (3*Compile_1090*Compile_11*Compile_5)/2.;
      std::complex<double> Compile_1099 = Compile_1045 + Compile_1046 + \
      Compile_1097 + Compile_1098;
      std::complex<double> Compile_1073 = Compile_1052 + Compile_16;
      std::complex<double> Compile_1074 = 1./Compile_1073;
      std::complex<double> Compile_1075 = mu3US*Compile_1074;
      std::complex<double> Compile_1076 = log(Compile_1075);
      std::complex<double> Compile_1077 = 0.5 + Compile_1076;
      std::complex<double> Compile_1039 = -0.25*(Compile_11*Compile_5);
      std::complex<double> Compile_1041 = Compile_1039 + Compile_1040;
      std::complex<double> Compile_1042 = pow(Compile_1041,2);
      std::complex<double> Compile_1091 = \
      -0.5*(Compile_1090*Compile_11*Compile_5);
      std::complex<double> Compile_1093 = \
      -0.5*(Compile_1019*Compile_5*Compile_6);
      std::complex<double> Compile_1094 = \
      -2*Compile_1090*Compile_11*Compile_5;
      std::complex<double> Compile_1095 = Compile_1093 + Compile_1094;
      std::complex<double> Compile_1096 = 3*Compile_1095;
      std::complex<double> Compile_1100 = -2*Compile_1099;
      std::complex<double> Compile_1101 = Compile_1096 + Compile_1100;
      std::complex<double> Compile_1102 = -0.015625*(Compile_1101*Compile_8);
      std::complex<double> Compile_1103 = \
      (-3*Compile_1077*Compile_1099*Compile_8)/16.;
      std::complex<double> Compile_1104 = Compile_1102 + Compile_1103;
      std::complex<double> Compile_1121 = \
      -0.0625*(Compile_117*Compile_68*Compile_8);
      std::complex<double> Compile_1122 = 1./sqrt(Compile_118);
      std::complex<double> Compile_1123 = (mu3US*Compile_1122)/2.;
      std::complex<double> Compile_1124 = log(Compile_1123);
      std::complex<double> Compile_1125 = 0.5 + Compile_1124;
      std::complex<double> Compile_1126 = \
      -0.125*(Compile_1125*Compile_117*Compile_68*Compile_8);
      std::complex<double> Compile_1127 = Compile_1121 + Compile_1126;
      std::complex<double> Compile_1128 = -2*Compile_1127;
      std::complex<double> Compile_1129 = Compile_1121 + Compile_1128;
      std::complex<double> Compile_1136 = pow(phi1,3);
      std::complex<double> Compile_1132 = pow(Compile_79,-2);
      std::complex<double> Compile_996 = pow(Compile_11,-2);
      std::complex<double> Compile_1134 = 2*m13dUS*phi1;
      std::complex<double> Compile_1135 = -2*m23dUS*phi1;
      std::complex<double> Compile_1137 = 2*lam1H3dUS*Compile_1136;
      std::complex<double> Compile_1138 = -(lam3H3dUS*Compile_1136);
      std::complex<double> Compile_1139 = -4*m12R3dUS*phi2;
      std::complex<double> Compile_1140 = -2*lam2H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1141 = lam3H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1142 = 2*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1143 = 2*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1144 = phi1*Compile_61;
      std::complex<double> Compile_1145 = Compile_1134 + Compile_1135 + \
      Compile_1137 + Compile_1138 + Compile_1139 + Compile_1140 + \
      Compile_1141 + Compile_1142 + Compile_1143 + Compile_1144;
      std::complex<double> Compile_1146 = pow(Compile_1145,2);
      std::complex<double> Compile_1149 = mu3US*Compile_1122;
      std::complex<double> Compile_1150 = log(Compile_1149);
      std::complex<double> Compile_1151 = 0.5 + Compile_1150;
      std::complex<double> Compile_1169 = -(Compile_117*Compile_68);
      std::complex<double> Compile_1155 = Compile_1051 + Compile_119;
      std::complex<double> Compile_1156 = 1./Compile_1155;
      std::complex<double> Compile_1157 = mu3US*Compile_1156;
      std::complex<double> Compile_1158 = log(Compile_1157);
      std::complex<double> Compile_1159 = 0.5 + Compile_1158;
      std::complex<double> Compile_1174 = Compile_1085 + Compile_1169;
      std::complex<double> Compile_1175 = pow(Compile_1174,2);
      std::complex<double> Compile_1165 = pow(g13dUS,4);
      std::complex<double> Compile_1166 = \
      (Compile_1019*Compile_5*Compile_6*Compile_8)/128.;
      std::complex<double> Compile_1167 = \
      (Compile_11*Compile_119*Compile_16*Compile_6*Compile_8)/128.;
      std::complex<double> Compile_1168 = \
      (Compile_11*Compile_119*Compile_18*Compile_5*Compile_8)/128.;
      std::complex<double> Compile_1170 = Compile_1040 + Compile_1085 + \
      Compile_1169;
      std::complex<double> Compile_1171 = \
      -0.015625*(Compile_1170*Compile_16*Compile_18*Compile_8);
      std::complex<double> Compile_1172 = pow(Compile_68,2);
      std::complex<double> Compile_1173 = \
      -0.0625*(Compile_1151*Compile_1172*Compile_120*Compile_8);
      std::complex<double> Compile_1176 = \
      (Compile_1159*Compile_1175*Compile_8)/16.;
      std::complex<double> Compile_1177 = Compile_1040 + Compile_1169;
      std::complex<double> Compile_1178 = pow(Compile_1177,2);
      std::complex<double> Compile_1179 = Compile_1052 + Compile_119;
      std::complex<double> Compile_1180 = 1./Compile_1179;
      std::complex<double> Compile_1181 = mu3US*Compile_1180;
      std::complex<double> Compile_1182 = log(Compile_1181);
      std::complex<double> Compile_1183 = 0.5 + Compile_1182;
      std::complex<double> Compile_1184 = \
      (Compile_1178*Compile_1183*Compile_8)/16.;
      std::complex<double> Compile_1187 = Compile_1051 + Compile_1052 + \
      Compile_119;
      std::complex<double> Compile_1188 = 1./Compile_1187;
      std::complex<double> Compile_1189 = mu3US*Compile_1188;
      std::complex<double> Compile_1190 = log(Compile_1189);
      std::complex<double> Compile_1191 = 0.5 + Compile_1190;
      std::complex<double> Compile_1133 = 1./Compile_11;
      std::complex<double> Compile_1201 = 1./pi;
      std::complex<double> Compile_1200 = \
      -0.0078125*(Compile_11*Compile_119*Compile_18*Compile_6*Compile_8);
      std::complex<double> Compile_1202 = \
      -0.03125*(Compile_11*Compile_1201*Compile_18*Compile_6);
      std::complex<double> Compile_1203 = \
      (Compile_11*Compile_119*Compile_1201*Compile_6)/16.;
      std::complex<double> Compile_1204 = Compile_1202 + Compile_1203;
      std::complex<double> Compile_1205 = \
      (Compile_119*Compile_1201*Compile_1204)/4.;
      std::complex<double> Compile_1206 = \
      -(Compile_11*Compile_117*Compile_6*Compile_68);
      std::complex<double> Compile_1207 = Compile_1046 + Compile_1206;
      std::complex<double> Compile_1208 = 2*Compile_119;
      std::complex<double> Compile_1209 = Compile_1052 + Compile_1208;
      std::complex<double> Compile_1210 = 1./Compile_1209;
      std::complex<double> Compile_1211 = mu3US*Compile_1210;
      std::complex<double> Compile_1212 = log(Compile_1211);
      std::complex<double> Compile_1213 = 0.5 + Compile_1212;
      std::complex<double> Compile_1214 = \
      -0.0625*(Compile_1207*Compile_1213*Compile_8);
      std::complex<double> Compile_1215 = Compile_1200 + Compile_1205 + \
      Compile_1214;
      std::complex<double> Compile_1221 = \
      -0.0625*(Compile_143*Compile_166*Compile_8);
      std::complex<double> Compile_1222 = 1./sqrt(Compile_167);
      std::complex<double> Compile_1223 = (mu3US*Compile_1222)/2.;
      std::complex<double> Compile_1224 = log(Compile_1223);
      std::complex<double> Compile_1225 = 0.5 + Compile_1224;
      std::complex<double> Compile_1226 = \
      -0.125*(Compile_1225*Compile_143*Compile_166*Compile_8);
      std::complex<double> Compile_1227 = Compile_1221 + Compile_1226;
      std::complex<double> Compile_1228 = -2*Compile_1227;
      std::complex<double> Compile_1229 = Compile_1221 + Compile_1228;
      std::complex<double> Compile_1147 = \
      (3*Compile_11*Compile_5*Compile_8)/128.;
      std::complex<double> Compile_1153 = (-3*Compile_11*Compile_5)/4.;
      std::complex<double> Compile_1232 = -2*m13dUS*phi1;
      std::complex<double> Compile_1233 = 2*m23dUS*phi1;
      std::complex<double> Compile_1234 = -2*lam1H3dUS*Compile_1136;
      std::complex<double> Compile_1235 = lam3H3dUS*Compile_1136;
      std::complex<double> Compile_1236 = 4*m12R3dUS*phi2;
      std::complex<double> Compile_1237 = 2*lam2H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1238 = -(lam3H3dUS*phi1*Compile_10);
      std::complex<double> Compile_1239 = -2*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1240 = -2*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1241 = Compile_1144 + Compile_1232 + \
      Compile_1233 + Compile_1234 + Compile_1235 + Compile_1236 + \
      Compile_1237 + Compile_1238 + Compile_1239 + Compile_1240;
      std::complex<double> Compile_1242 = pow(Compile_1241,2);
      std::complex<double> Compile_1244 = mu3US*Compile_1222;
      std::complex<double> Compile_1245 = log(Compile_1244);
      std::complex<double> Compile_1246 = 0.5 + Compile_1245;
      std::complex<double> Compile_1261 = -(Compile_143*Compile_166);
      std::complex<double> Compile_1249 = Compile_1051 + Compile_168;
      std::complex<double> Compile_1250 = 1./Compile_1249;
      std::complex<double> Compile_1251 = mu3US*Compile_1250;
      std::complex<double> Compile_1252 = log(Compile_1251);
      std::complex<double> Compile_1253 = 0.5 + Compile_1252;
      std::complex<double> Compile_1266 = Compile_1085 + Compile_1261;
      std::complex<double> Compile_1267 = pow(Compile_1266,2);
      std::complex<double> Compile_1259 = \
      (Compile_11*Compile_16*Compile_168*Compile_6*Compile_8)/128.;
      std::complex<double> Compile_1260 = \
      (Compile_11*Compile_168*Compile_18*Compile_5*Compile_8)/128.;
      std::complex<double> Compile_1262 = Compile_1040 + Compile_1085 + \
      Compile_1261;
      std::complex<double> Compile_1263 = \
      -0.015625*(Compile_1262*Compile_16*Compile_18*Compile_8);
      std::complex<double> Compile_1264 = pow(Compile_143,2);
      std::complex<double> Compile_1265 = \
      -0.0625*(Compile_1246*Compile_1264*Compile_169*Compile_8);
      std::complex<double> Compile_1268 = \
      (Compile_1253*Compile_1267*Compile_8)/16.;
      std::complex<double> Compile_1269 = Compile_1040 + Compile_1261;
      std::complex<double> Compile_1270 = pow(Compile_1269,2);
      std::complex<double> Compile_1271 = Compile_1052 + Compile_168;
      std::complex<double> Compile_1272 = 1./Compile_1271;
      std::complex<double> Compile_1273 = mu3US*Compile_1272;
      std::complex<double> Compile_1274 = log(Compile_1273);
      std::complex<double> Compile_1275 = 0.5 + Compile_1274;
      std::complex<double> Compile_1276 = \
      (Compile_1270*Compile_1275*Compile_8)/16.;
      std::complex<double> Compile_1279 = Compile_1051 + Compile_1052 + \
      Compile_168;
      std::complex<double> Compile_1280 = 1./Compile_1279;
      std::complex<double> Compile_1281 = mu3US*Compile_1280;
      std::complex<double> Compile_1282 = log(Compile_1281);
      std::complex<double> Compile_1283 = 0.5 + Compile_1282;
      std::complex<double> Compile_1292 = \
      -0.0078125*(Compile_11*Compile_168*Compile_18*Compile_6*Compile_8);
      std::complex<double> Compile_1293 = \
      (Compile_11*Compile_1201*Compile_168*Compile_6)/16.;
      std::complex<double> Compile_1294 = Compile_1202 + Compile_1293;
      std::complex<double> Compile_1295 = \
      (Compile_1201*Compile_1294*Compile_168)/4.;
      std::complex<double> Compile_1296 = \
      -(Compile_11*Compile_143*Compile_166*Compile_6);
      std::complex<double> Compile_1297 = Compile_1046 + Compile_1296;
      std::complex<double> Compile_1298 = 2*Compile_168;
      std::complex<double> Compile_1299 = Compile_1052 + Compile_1298;
      std::complex<double> Compile_1300 = 1./Compile_1299;
      std::complex<double> Compile_1301 = mu3US*Compile_1300;
      std::complex<double> Compile_1302 = log(Compile_1301);
      std::complex<double> Compile_1303 = 0.5 + Compile_1302;
      std::complex<double> Compile_1304 = \
      -0.0625*(Compile_1297*Compile_1303*Compile_8);
      std::complex<double> Compile_1305 = Compile_1292 + Compile_1295 + \
      Compile_1304;
      std::complex<double> Compile_1311 = -(Compile_330*Compile_348);
      std::complex<double> Compile_1316 = Compile_118 + Compile_1311;
      std::complex<double> Compile_1317 = pow(Compile_1316,2);
      std::complex<double> Compile_1308 = Compile_1085 + Compile_1169 + \
      Compile_362;
      std::complex<double> Compile_1309 = \
      -0.03125*(Compile_1308*Compile_16*Compile_363*Compile_8);
      std::complex<double> Compile_1310 = \
      (Compile_11*Compile_1201*Compile_363*Compile_5)/16.;
      std::complex<double> Compile_1312 = Compile_1085 + Compile_118 + \
      Compile_1311;
      std::complex<double> Compile_1313 = \
      -0.125*(Compile_1201*Compile_1312*Compile_16);
      std::complex<double> Compile_1314 = Compile_1310 + Compile_1313;
      std::complex<double> Compile_1315 = \
      (Compile_119*Compile_1201*Compile_1314)/4.;
      std::complex<double> Compile_1318 = Compile_119 + Compile_363;
      std::complex<double> Compile_1319 = 1./Compile_1318;
      std::complex<double> Compile_1320 = mu3US*Compile_1319;
      std::complex<double> Compile_1321 = log(Compile_1320);
      std::complex<double> Compile_1322 = 0.5 + Compile_1321;
      std::complex<double> Compile_1323 = \
      (Compile_1317*Compile_1322*Compile_8)/16.;
      std::complex<double> Compile_1324 = Compile_118 + Compile_362;
      std::complex<double> Compile_1325 = \
      -0.5*(Compile_11*Compile_1324*Compile_5);
      std::complex<double> Compile_1326 = Compile_1044 + Compile_1317 + \
      Compile_1325;
      std::complex<double> Compile_1327 = Compile_1051 + Compile_119 + \
      Compile_363;
      std::complex<double> Compile_1328 = 1./Compile_1327;
      std::complex<double> Compile_1329 = mu3US*Compile_1328;
      std::complex<double> Compile_1330 = log(Compile_1329);
      std::complex<double> Compile_1331 = 0.5 + Compile_1330;
      std::complex<double> Compile_1332 = \
      -0.0625*(Compile_1326*Compile_1331*Compile_8);
      std::complex<double> Compile_1333 = Compile_1309 + Compile_1315 + \
      Compile_1323 + Compile_1332;
      std::complex<double> Compile_1341 = Compile_1169 + Compile_362;
      std::complex<double> Compile_1342 = pow(Compile_1341,2);
      std::complex<double> Compile_1336 = \
      -0.03125*(Compile_119*Compile_1312*Compile_16*Compile_8);
      std::complex<double> Compile_1337 = \
      (Compile_11*Compile_119*Compile_1201*Compile_5)/16.;
      std::complex<double> Compile_1338 = \
      -0.125*(Compile_1201*Compile_1308*Compile_16);
      std::complex<double> Compile_1339 = Compile_1337 + Compile_1338;
      std::complex<double> Compile_1340 = \
      (Compile_1201*Compile_1339*Compile_363)/4.;
      std::complex<double> Compile_1343 = \
      (Compile_1322*Compile_1342*Compile_8)/16.;
      std::complex<double> Compile_1344 = Compile_1044 + Compile_1325 + \
      Compile_1342;
      std::complex<double> Compile_1345 = \
      -0.0625*(Compile_1331*Compile_1344*Compile_8);
      std::complex<double> Compile_1346 = Compile_1336 + Compile_1340 + \
      Compile_1343 + Compile_1345;
      std::complex<double> Compile_1355 = Compile_1311 + Compile_167;
      std::complex<double> Compile_1356 = pow(Compile_1355,2);
      std::complex<double> Compile_1349 = Compile_1085 + Compile_1261 + \
      Compile_362;
      std::complex<double> Compile_1350 = \
      -0.03125*(Compile_1349*Compile_16*Compile_363*Compile_8);
      std::complex<double> Compile_1351 = Compile_1085 + Compile_1311 + \
      Compile_167;
      std::complex<double> Compile_1352 = \
      -0.125*(Compile_1201*Compile_1351*Compile_16);
      std::complex<double> Compile_1353 = Compile_1310 + Compile_1352;
      std::complex<double> Compile_1354 = \
      (Compile_1201*Compile_1353*Compile_168)/4.;
      std::complex<double> Compile_1357 = Compile_168 + Compile_363;
      std::complex<double> Compile_1358 = 1./Compile_1357;
      std::complex<double> Compile_1359 = mu3US*Compile_1358;
      std::complex<double> Compile_1360 = log(Compile_1359);
      std::complex<double> Compile_1361 = 0.5 + Compile_1360;
      std::complex<double> Compile_1362 = \
      (Compile_1356*Compile_1361*Compile_8)/16.;
      std::complex<double> Compile_1363 = Compile_167 + Compile_362;
      std::complex<double> Compile_1364 = \
      -0.5*(Compile_11*Compile_1363*Compile_5);
      std::complex<double> Compile_1365 = Compile_1044 + Compile_1356 + \
      Compile_1364;
      std::complex<double> Compile_1366 = Compile_1051 + Compile_168 + \
      Compile_363;
      std::complex<double> Compile_1367 = 1./Compile_1366;
      std::complex<double> Compile_1368 = mu3US*Compile_1367;
      std::complex<double> Compile_1369 = log(Compile_1368);
      std::complex<double> Compile_1370 = 0.5 + Compile_1369;
      std::complex<double> Compile_1371 = \
      -0.0625*(Compile_1365*Compile_1370*Compile_8);
      std::complex<double> Compile_1372 = Compile_1350 + Compile_1354 + \
      Compile_1362 + Compile_1371;
      std::complex<double> Compile_1380 = Compile_1261 + Compile_362;
      std::complex<double> Compile_1381 = pow(Compile_1380,2);
      std::complex<double> Compile_1375 = \
      -0.03125*(Compile_1351*Compile_16*Compile_168*Compile_8);
      std::complex<double> Compile_1376 = \
      (Compile_11*Compile_1201*Compile_168*Compile_5)/16.;
      std::complex<double> Compile_1377 = \
      -0.125*(Compile_1201*Compile_1349*Compile_16);
      std::complex<double> Compile_1378 = Compile_1376 + Compile_1377;
      std::complex<double> Compile_1379 = \
      (Compile_1201*Compile_1378*Compile_363)/4.;
      std::complex<double> Compile_1382 = \
      (Compile_1361*Compile_1381*Compile_8)/16.;
      std::complex<double> Compile_1383 = Compile_1044 + Compile_1364 + \
      Compile_1381;
      std::complex<double> Compile_1384 = \
      -0.0625*(Compile_1370*Compile_1383*Compile_8);
      std::complex<double> Compile_1385 = Compile_1375 + Compile_1379 + \
      Compile_1382 + Compile_1384;
      std::complex<double> Compile_1397 = pow(phi2,3);
      std::complex<double> Compile_1414 = -(Compile_435*Compile_454);
      std::complex<double> Compile_1419 = Compile_118 + Compile_1414;
      std::complex<double> Compile_1420 = pow(Compile_1419,2);
      std::complex<double> Compile_1411 = Compile_1085 + Compile_1169 + \
      Compile_468;
      std::complex<double> Compile_1412 = \
      -0.03125*(Compile_1411*Compile_16*Compile_469*Compile_8);
      std::complex<double> Compile_1413 = \
      (Compile_11*Compile_1201*Compile_469*Compile_5)/16.;
      std::complex<double> Compile_1415 = Compile_1085 + Compile_118 + \
      Compile_1414;
      std::complex<double> Compile_1416 = \
      -0.125*(Compile_1201*Compile_1415*Compile_16);
      std::complex<double> Compile_1417 = Compile_1413 + Compile_1416;
      std::complex<double> Compile_1418 = \
      (Compile_119*Compile_1201*Compile_1417)/4.;
      std::complex<double> Compile_1421 = Compile_119 + Compile_469;
      std::complex<double> Compile_1422 = 1./Compile_1421;
      std::complex<double> Compile_1423 = mu3US*Compile_1422;
      std::complex<double> Compile_1424 = log(Compile_1423);
      std::complex<double> Compile_1425 = 0.5 + Compile_1424;
      std::complex<double> Compile_1426 = \
      (Compile_1420*Compile_1425*Compile_8)/16.;
      std::complex<double> Compile_1427 = Compile_118 + Compile_468;
      std::complex<double> Compile_1428 = \
      -0.5*(Compile_11*Compile_1427*Compile_5);
      std::complex<double> Compile_1429 = Compile_1044 + Compile_1420 + \
      Compile_1428;
      std::complex<double> Compile_1430 = Compile_1051 + Compile_119 + \
      Compile_469;
      std::complex<double> Compile_1431 = 1./Compile_1430;
      std::complex<double> Compile_1432 = mu3US*Compile_1431;
      std::complex<double> Compile_1433 = log(Compile_1432);
      std::complex<double> Compile_1434 = 0.5 + Compile_1433;
      std::complex<double> Compile_1435 = \
      -0.0625*(Compile_1429*Compile_1434*Compile_8);
      std::complex<double> Compile_1436 = Compile_1412 + Compile_1418 + \
      Compile_1426 + Compile_1435;
      std::complex<double> Compile_1443 = Compile_1169 + Compile_468;
      std::complex<double> Compile_1444 = pow(Compile_1443,2);
      std::complex<double> Compile_1439 = \
      -0.03125*(Compile_119*Compile_1415*Compile_16*Compile_8);
      std::complex<double> Compile_1440 = \
      -0.125*(Compile_1201*Compile_1411*Compile_16);
      std::complex<double> Compile_1441 = Compile_1337 + Compile_1440;
      std::complex<double> Compile_1442 = \
      (Compile_1201*Compile_1441*Compile_469)/4.;
      std::complex<double> Compile_1445 = \
      (Compile_1425*Compile_1444*Compile_8)/16.;
      std::complex<double> Compile_1446 = Compile_1044 + Compile_1428 + \
      Compile_1444;
      std::complex<double> Compile_1447 = \
      -0.0625*(Compile_1434*Compile_1446*Compile_8);
      std::complex<double> Compile_1448 = Compile_1439 + Compile_1442 + \
      Compile_1445 + Compile_1447;
      std::complex<double> Compile_1457 = Compile_1414 + Compile_167;
      std::complex<double> Compile_1458 = pow(Compile_1457,2);
      std::complex<double> Compile_1451 = Compile_1085 + Compile_1261 + \
      Compile_468;
      std::complex<double> Compile_1452 = \
      -0.03125*(Compile_1451*Compile_16*Compile_469*Compile_8);
      std::complex<double> Compile_1453 = Compile_1085 + Compile_1414 + \
      Compile_167;
      std::complex<double> Compile_1454 = \
      -0.125*(Compile_1201*Compile_1453*Compile_16);
      std::complex<double> Compile_1455 = Compile_1413 + Compile_1454;
      std::complex<double> Compile_1456 = \
      (Compile_1201*Compile_1455*Compile_168)/4.;
      std::complex<double> Compile_1459 = Compile_168 + Compile_469;
      std::complex<double> Compile_1460 = 1./Compile_1459;
      std::complex<double> Compile_1461 = mu3US*Compile_1460;
      std::complex<double> Compile_1462 = log(Compile_1461);
      std::complex<double> Compile_1463 = 0.5 + Compile_1462;
      std::complex<double> Compile_1464 = \
      (Compile_1458*Compile_1463*Compile_8)/16.;
      std::complex<double> Compile_1465 = Compile_167 + Compile_468;
      std::complex<double> Compile_1466 = \
      -0.5*(Compile_11*Compile_1465*Compile_5);
      std::complex<double> Compile_1467 = Compile_1044 + Compile_1458 + \
      Compile_1466;
      std::complex<double> Compile_1468 = Compile_1051 + Compile_168 + \
      Compile_469;
      std::complex<double> Compile_1469 = 1./Compile_1468;
      std::complex<double> Compile_1470 = mu3US*Compile_1469;
      std::complex<double> Compile_1471 = log(Compile_1470);
      std::complex<double> Compile_1472 = 0.5 + Compile_1471;
      std::complex<double> Compile_1473 = \
      -0.0625*(Compile_1467*Compile_1472*Compile_8);
      std::complex<double> Compile_1474 = Compile_1452 + Compile_1456 + \
      Compile_1464 + Compile_1473;
      std::complex<double> Compile_1481 = Compile_1261 + Compile_468;
      std::complex<double> Compile_1482 = pow(Compile_1481,2);
      std::complex<double> Compile_1477 = \
      -0.03125*(Compile_1453*Compile_16*Compile_168*Compile_8);
      std::complex<double> Compile_1478 = \
      -0.125*(Compile_1201*Compile_1451*Compile_16);
      std::complex<double> Compile_1479 = Compile_1376 + Compile_1478;
      std::complex<double> Compile_1480 = \
      (Compile_1201*Compile_1479*Compile_469)/4.;
      std::complex<double> Compile_1483 = \
      (Compile_1463*Compile_1482*Compile_8)/16.;
      std::complex<double> Compile_1484 = Compile_1044 + Compile_1466 + \
      Compile_1482;
      std::complex<double> Compile_1485 = \
      -0.0625*(Compile_1472*Compile_1484*Compile_8);
      std::complex<double> Compile_1486 = Compile_1477 + Compile_1480 + \
      Compile_1483 + Compile_1485;
      std::complex<double> Compile_1388 = lam4H3dUS + Compile_179;
      std::complex<double> Compile_1389 = pow(Compile_1388,2);
      std::complex<double> Compile_1390 = 4*m12R3dUS*phi1;
      std::complex<double> Compile_1391 = 2*m13dUS*phi2;
      std::complex<double> Compile_1392 = -2*m23dUS*phi2;
      std::complex<double> Compile_1393 = 2*lam1H3dUS*phi2*Compile_9;
      std::complex<double> Compile_1394 = -(lam3H3dUS*phi2*Compile_9);
      std::complex<double> Compile_1395 = -(lam4H3dUS*phi2*Compile_9);
      std::complex<double> Compile_1396 = -3*lam5H3dUS*phi2*Compile_9;
      std::complex<double> Compile_1398 = -2*lam2H3dUS*Compile_1397;
      std::complex<double> Compile_1399 = lam3H3dUS*Compile_1397;
      std::complex<double> Compile_1400 = lam4H3dUS*Compile_1397;
      std::complex<double> Compile_1401 = -(lam5H3dUS*Compile_1397);
      std::complex<double> Compile_1499 = -2*m13dUS*phi2;
      std::complex<double> Compile_1500 = 2*m23dUS*phi2;
      std::complex<double> Compile_1501 = -6*lam1H3dUS*phi2*Compile_9;
      std::complex<double> Compile_1502 = -4*phi2*Compile_71*Compile_9;
      std::complex<double> Compile_1503 = 6*lam2H3dUS*Compile_1397;
      std::complex<double> Compile_1504 = \
      phi2*Compile_137*Compile_138*Compile_71;
      std::complex<double> Compile_1505 = phi2*Compile_608;
      std::complex<double> Compile_1506 = Compile_1390 + Compile_1499 + \
      Compile_1500 + Compile_1501 + Compile_1502 + Compile_1503 + \
      Compile_1504 + Compile_1505;
      std::complex<double> Compile_632 = pow(Compile_631,-3);
      std::complex<double> Compile_1523 = 1./sqrt(Compile_654);
      std::complex<double> Compile_1541 = -(Compile_619*Compile_637);
      std::complex<double> Compile_1550 = Compile_1085 + Compile_1541;
      std::complex<double> Compile_1551 = pow(Compile_1550,2);
      std::complex<double> Compile_1528 = -6*lam1H3dUS*Compile_1136;
      std::complex<double> Compile_1529 = lam4H3dUS*Compile_1136;
      std::complex<double> Compile_1530 = lam5H3dUS*Compile_1136;
      std::complex<double> Compile_1531 = 6*lam2H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1532 = -5*lam3H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1533 = -5*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1534 = -5*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_1535 = phi1*Compile_608;
      std::complex<double> Compile_1536 = Compile_1232 + Compile_1233 + \
      Compile_1235 + Compile_1236 + Compile_1528 + Compile_1529 + \
      Compile_1530 + Compile_1531 + Compile_1532 + Compile_1533 + \
      Compile_1534 + Compile_1535;
      std::complex<double> Compile_1537 = pow(Compile_1536,2);
      std::complex<double> Compile_1544 = pow(Compile_619,2);
      std::complex<double> Compile_1545 = pow(Compile_631,-2);
      std::complex<double> Compile_1546 = mu3US*Compile_1523;
      std::complex<double> Compile_1547 = log(Compile_1546);
      std::complex<double> Compile_1548 = 0.5 + Compile_1547;
      std::complex<double> Compile_1549 = \
      -0.0625*(Compile_1544*Compile_1545*Compile_1548*Compile_8);
      std::complex<double> Compile_1574 = Compile_1040 + Compile_1541;
      std::complex<double> Compile_1575 = pow(Compile_1574,2);
      std::complex<double> Compile_634 = \
      -(Compile_344*Compile_616*Compile_627*Compile_80);
      std::complex<double> Compile_635 = -8 + Compile_634;
      std::complex<double> Compile_636 = pow(Compile_635,2);
      std::complex<double> Compile_1600 = Compile_118 + Compile_1541;
      std::complex<double> Compile_1601 = pow(Compile_1600,2);
      std::complex<double> Compile_651 = \
      Compile_344*Compile_616*Compile_627*Compile_80;
      std::complex<double> Compile_652 = 8 + Compile_651;
      std::complex<double> Compile_653 = pow(Compile_652,2);
      std::complex<double> Compile_1596 = Compile_1085 + Compile_118 + \
      Compile_1541;
      std::complex<double> Compile_1593 = Compile_1085 + Compile_1169 + \
      Compile_654;
      std::complex<double> Compile_1602 = Compile_119 + Compile_655;
      std::complex<double> Compile_1603 = 1./Compile_1602;
      std::complex<double> Compile_1604 = mu3US*Compile_1603;
      std::complex<double> Compile_1605 = log(Compile_1604);
      std::complex<double> Compile_1606 = 0.5 + Compile_1605;
      std::complex<double> Compile_1623 = Compile_1169 + Compile_654;
      std::complex<double> Compile_1624 = pow(Compile_1623,2);
      std::complex<double> Compile_1608 = Compile_118 + Compile_654;
      std::complex<double> Compile_1609 = \
      -0.5*(Compile_11*Compile_1608*Compile_5);
      std::complex<double> Compile_1611 = Compile_1051 + Compile_119 + \
      Compile_655;
      std::complex<double> Compile_1612 = 1./Compile_1611;
      std::complex<double> Compile_1613 = mu3US*Compile_1612;
      std::complex<double> Compile_1614 = log(Compile_1613);
      std::complex<double> Compile_1615 = 0.5 + Compile_1614;
      std::complex<double> Compile_1640 = phi2*Compile_21*Compile_79;
      std::complex<double> Compile_1641 = -2*lam1H3dUS*phi1*Compile_62;
      std::complex<double> Compile_1642 = Compile_1640 + Compile_1641;
      std::complex<double> Compile_1634 = phi1*Compile_21*Compile_79;
      std::complex<double> Compile_1635 = -(lam3H3dUS*phi2*Compile_62);
      std::complex<double> Compile_1636 = Compile_1634 + Compile_1635;
      std::complex<double> Compile_1630 = 8*lam2H3dUS*phi2*Compile_79;
      std::complex<double> Compile_1631 = -(phi1*Compile_21*Compile_62);
      std::complex<double> Compile_1632 = Compile_1630 + Compile_1631;
      std::complex<double> Compile_1644 = 4*lam3H3dUS*phi1*Compile_79;
      std::complex<double> Compile_1645 = -(phi2*Compile_21*Compile_62);
      std::complex<double> Compile_1646 = Compile_1644 + Compile_1645;
      std::complex<double> Compile_1652 = Compile_1208 + Compile_655;
      std::complex<double> Compile_1653 = 1./Compile_1652;
      std::complex<double> Compile_1654 = mu3US*Compile_1653;
      std::complex<double> Compile_1655 = log(Compile_1654);
      std::complex<double> Compile_1656 = 0.5 + Compile_1655;
      std::complex<double> Compile_639 = \
      -(Compile_140*Compile_616*Compile_627*Compile_80);
      std::complex<double> Compile_640 = -8 + Compile_639;
      std::complex<double> Compile_641 = pow(Compile_640,2);
      std::complex<double> Compile_1595 = \
      (Compile_11*Compile_1201*Compile_5*Compile_655)/16.;
      std::complex<double> Compile_1692 = Compile_1541 + Compile_167;
      std::complex<double> Compile_1693 = pow(Compile_1692,2);
      std::complex<double> Compile_694 = \
      Compile_140*Compile_616*Compile_627*Compile_80;
      std::complex<double> Compile_695 = 8 + Compile_694;
      std::complex<double> Compile_696 = pow(Compile_695,2);
      std::complex<double> Compile_1688 = Compile_1085 + Compile_1541 + \
      Compile_167;
      std::complex<double> Compile_1686 = Compile_1085 + Compile_1261 + \
      Compile_654;
      std::complex<double> Compile_1694 = Compile_168 + Compile_655;
      std::complex<double> Compile_1695 = 1./Compile_1694;
      std::complex<double> Compile_1696 = mu3US*Compile_1695;
      std::complex<double> Compile_1697 = log(Compile_1696);
      std::complex<double> Compile_1698 = 0.5 + Compile_1697;
      std::complex<double> Compile_1715 = Compile_1261 + Compile_654;
      std::complex<double> Compile_1716 = pow(Compile_1715,2);
      std::complex<double> Compile_1700 = Compile_167 + Compile_654;
      std::complex<double> Compile_1701 = \
      -0.5*(Compile_11*Compile_1700*Compile_5);
      std::complex<double> Compile_1703 = Compile_1051 + Compile_168 + \
      Compile_655;
      std::complex<double> Compile_1704 = 1./Compile_1703;
      std::complex<double> Compile_1705 = mu3US*Compile_1704;
      std::complex<double> Compile_1706 = log(Compile_1705);
      std::complex<double> Compile_1707 = 0.5 + Compile_1706;
      std::complex<double> Compile_1810 = pow(phi2,5);
      std::complex<double> Compile_1853 = -2*lam1H3dUS*phi1*Compile_154;
      std::complex<double> Compile_1854 = Compile_1640 + Compile_1853;
      std::complex<double> Compile_1848 = -(lam3H3dUS*phi2*Compile_154);
      std::complex<double> Compile_1849 = Compile_1634 + Compile_1848;
      std::complex<double> Compile_1845 = -(phi1*Compile_154*Compile_21);
      std::complex<double> Compile_1846 = Compile_1630 + Compile_1845;
      std::complex<double> Compile_1856 = phi2*Compile_126*Compile_21;
      std::complex<double> Compile_1857 = Compile_1644 + Compile_1856;
      std::complex<double> Compile_1863 = Compile_1298 + Compile_655;
      std::complex<double> Compile_1864 = 1./Compile_1863;
      std::complex<double> Compile_1865 = mu3US*Compile_1864;
      std::complex<double> Compile_1866 = log(Compile_1865);
      std::complex<double> Compile_1867 = 0.5 + Compile_1866;
      std::complex<double> Compile_1669 = \
      -(Compile_1506*Compile_21*Compile_79);
      std::complex<double> Compile_1670 = 2*lam3H3dUS*phi2*Compile_577;
      std::complex<double> Compile_1671 = -(lam1H3dUS*phi1*Compile_661);
      std::complex<double> Compile_1672 = Compile_1670 + Compile_1671;
      std::complex<double> Compile_1677 = 8*lam2H3dUS*phi2*Compile_577;
      std::complex<double> Compile_1678 = -(lam3H3dUS*phi1*Compile_661);
      std::complex<double> Compile_1679 = Compile_1677 + Compile_1678;
      std::complex<double> Compile_1680 = 4*Compile_1679*Compile_79;
      std::complex<double> Compile_1896 = Compile_1541 + Compile_362;
      std::complex<double> Compile_1897 = pow(Compile_1896,2);
      std::complex<double> Compile_1892 = Compile_1040 + Compile_1541 + \
      Compile_362;
      std::complex<double> Compile_1889 = Compile_1040 + Compile_1311 + \
      Compile_654;
      std::complex<double> Compile_1898 = Compile_363 + Compile_655;
      std::complex<double> Compile_1899 = 1./Compile_1898;
      std::complex<double> Compile_1900 = mu3US*Compile_1899;
      std::complex<double> Compile_1901 = log(Compile_1900);
      std::complex<double> Compile_1902 = 0.5 + Compile_1901;
      std::complex<double> Compile_1920 = Compile_1311 + Compile_654;
      std::complex<double> Compile_1921 = pow(Compile_1920,2);
      std::complex<double> Compile_1904 = Compile_362 + Compile_654;
      std::complex<double> Compile_1905 = \
      -0.5*(Compile_11*Compile_1904*Compile_6);
      std::complex<double> Compile_1907 = Compile_1052 + Compile_363 + \
      Compile_655;
      std::complex<double> Compile_1908 = 1./Compile_1907;
      std::complex<double> Compile_1909 = mu3US*Compile_1908;
      std::complex<double> Compile_1910 = log(Compile_1909);
      std::complex<double> Compile_1911 = 0.5 + Compile_1910;
      std::complex<double> Compile_1937 = 2*lam5H3dUS*phi2*Compile_291;
      std::complex<double> Compile_1938 = -(lam1H3dUS*phi1*Compile_369);
      std::complex<double> Compile_1939 = Compile_1937 + Compile_1938;
      std::complex<double> Compile_1931 = 4*lam5H3dUS*phi1*Compile_291;
      std::complex<double> Compile_1932 = -(phi2*Compile_180*Compile_369);
      std::complex<double> Compile_1933 = Compile_1931 + Compile_1932;
      std::complex<double> Compile_1927 = 8*lam2H3dUS*phi2*Compile_291;
      std::complex<double> Compile_1928 = -(lam5H3dUS*phi1*Compile_369);
      std::complex<double> Compile_1929 = Compile_1927 + Compile_1928;
      std::complex<double> Compile_1941 = 4*phi1*Compile_180*Compile_291;
      std::complex<double> Compile_1942 = -(lam5H3dUS*phi2*Compile_369);
      std::complex<double> Compile_1943 = Compile_1941 + Compile_1942;
      std::complex<double> Compile_1949 = pow(Compile_341,-2);
      std::complex<double> Compile_1950 = 2*Compile_363;
      std::complex<double> Compile_1951 = Compile_1950 + Compile_655;
      std::complex<double> Compile_1952 = 1./Compile_1951;
      std::complex<double> Compile_1953 = mu3US*Compile_1952;
      std::complex<double> Compile_1954 = log(Compile_1953);
      std::complex<double> Compile_1955 = 0.5 + Compile_1954;
      std::complex<double> Compile_1891 = \
      (Compile_11*Compile_1201*Compile_6*Compile_655)/16.;
      std::complex<double> Compile_1989 = Compile_1541 + Compile_468;
      std::complex<double> Compile_1990 = pow(Compile_1989,2);
      std::complex<double> Compile_1985 = Compile_1040 + Compile_1541 + \
      Compile_468;
      std::complex<double> Compile_1983 = Compile_1040 + Compile_1414 + \
      Compile_654;
      std::complex<double> Compile_1991 = Compile_469 + Compile_655;
      std::complex<double> Compile_1992 = 1./Compile_1991;
      std::complex<double> Compile_1993 = mu3US*Compile_1992;
      std::complex<double> Compile_1994 = log(Compile_1993);
      std::complex<double> Compile_1995 = 0.5 + Compile_1994;
      std::complex<double> Compile_2013 = Compile_1414 + Compile_654;
      std::complex<double> Compile_2014 = pow(Compile_2013,2);
      std::complex<double> Compile_1997 = Compile_468 + Compile_654;
      std::complex<double> Compile_1998 = \
      -0.5*(Compile_11*Compile_1997*Compile_6);
      std::complex<double> Compile_2000 = Compile_1052 + Compile_469 + \
      Compile_655;
      std::complex<double> Compile_2001 = 1./Compile_2000;
      std::complex<double> Compile_2002 = mu3US*Compile_2001;
      std::complex<double> Compile_2003 = log(Compile_2002);
      std::complex<double> Compile_2004 = 0.5 + Compile_2003;
      std::complex<double> Compile_1804 = \
      4*lam1H3dUS*lam4H3dUS*lam5H3dUS*Compile_1397*Compile_9;
      std::complex<double> Compile_2126 = -(lam1H3dUS*phi1*Compile_431);
      std::complex<double> Compile_2127 = Compile_1937 + Compile_2126;
      std::complex<double> Compile_2121 = -(phi2*Compile_180*Compile_431);
      std::complex<double> Compile_2122 = Compile_1931 + Compile_2121;
      std::complex<double> Compile_2118 = -(lam5H3dUS*phi1*Compile_431);
      std::complex<double> Compile_2119 = Compile_1927 + Compile_2118;
      std::complex<double> Compile_2129 = -(lam5H3dUS*phi2*Compile_431);
      std::complex<double> Compile_2130 = Compile_1941 + Compile_2129;
      std::complex<double> Compile_2136 = pow(Compile_447,-2);
      std::complex<double> Compile_2137 = 2*Compile_469;
      std::complex<double> Compile_2138 = Compile_2137 + Compile_655;
      std::complex<double> Compile_2139 = 1./Compile_2138;
      std::complex<double> Compile_2140 = mu3US*Compile_2139;
      std::complex<double> Compile_2141 = log(Compile_2140);
      std::complex<double> Compile_2142 = 0.5 + Compile_2141;
      std::complex<double> Compile_1968 = \
      -2*lam5H3dUS*Compile_1506*Compile_291;
      std::complex<double> Compile_1969 = 2*phi2*Compile_180*Compile_577;
      std::complex<double> Compile_1970 = Compile_1671 + Compile_1969;
      std::complex<double> Compile_1975 = -(phi1*Compile_180*Compile_661);
      std::complex<double> Compile_1976 = Compile_1677 + Compile_1975;
      std::complex<double> Compile_1977 = 4*Compile_1976*Compile_291;
      std::complex<double> Compile_1498 = pow(Compile_577,-6);
      std::complex<double> Compile_1508 = 2*phi2*Compile_577*Compile_71;
      std::complex<double> Compile_1515 = 24*lam2H3dUS*phi2*Compile_577;
      std::complex<double> Compile_2167 = 4*phi1*Compile_577;
      std::complex<double> Compile_2168 = -(phi2*Compile_788);
      std::complex<double> Compile_2169 = Compile_2167 + Compile_2168;
      std::complex<double> Compile_803 = pow(Compile_802,-3);
      std::complex<double> Compile_1538 = \
      (Compile_1019*Compile_3*Compile_8)/128.;
      std::complex<double> Compile_1540 = (Compile_11*Compile_5)/2.;
      std::complex<double> Compile_2181 = 1./sqrt(Compile_826);
      std::complex<double> Compile_2190 = -(Compile_792*Compile_809);
      std::complex<double> Compile_1558 = (7*Compile_1019*Compile_3)/16.;
      std::complex<double> Compile_2199 = Compile_1085 + Compile_2190;
      std::complex<double> Compile_2200 = pow(Compile_2199,2);
      std::complex<double> Compile_2186 = -(phi1*Compile_608);
      std::complex<double> Compile_2187 = Compile_1232 + Compile_1233 + \
      Compile_1235 + Compile_1236 + Compile_1528 + Compile_1529 + \
      Compile_1530 + Compile_1531 + Compile_1532 + Compile_1533 + \
      Compile_1534 + Compile_2186;
      std::complex<double> Compile_2188 = pow(Compile_2187,2);
      std::complex<double> Compile_1569 = \
      (Compile_1015*Compile_1019*Compile_8)/128.;
      std::complex<double> Compile_1571 = (Compile_11*Compile_6)/2.;
      std::complex<double> Compile_2193 = pow(Compile_792,2);
      std::complex<double> Compile_2194 = pow(Compile_802,-2);
      std::complex<double> Compile_2195 = mu3US*Compile_2181;
      std::complex<double> Compile_2196 = log(Compile_2195);
      std::complex<double> Compile_2197 = 0.5 + Compile_2196;
      std::complex<double> Compile_2198 = \
      -0.0625*(Compile_2193*Compile_2194*Compile_2197*Compile_8);
      std::complex<double> Compile_1582 = \
      (7*Compile_1015*Compile_1019)/16.;
      std::complex<double> Compile_2220 = Compile_1040 + Compile_2190;
      std::complex<double> Compile_2221 = pow(Compile_2220,2);
      std::complex<double> Compile_806 = \
      -(Compile_344*Compile_627*Compile_80*Compile_805);
      std::complex<double> Compile_807 = -8 + Compile_806;
      std::complex<double> Compile_808 = pow(Compile_807,2);
      std::complex<double> Compile_2245 = Compile_118 + Compile_2190;
      std::complex<double> Compile_2246 = pow(Compile_2245,2);
      std::complex<double> Compile_823 = \
      Compile_344*Compile_627*Compile_80*Compile_805;
      std::complex<double> Compile_824 = 8 + Compile_823;
      std::complex<double> Compile_825 = pow(Compile_824,2);
      std::complex<double> Compile_2241 = Compile_1085 + Compile_118 + \
      Compile_2190;
      std::complex<double> Compile_2238 = Compile_1085 + Compile_1169 + \
      Compile_826;
      std::complex<double> Compile_2247 = Compile_119 + Compile_827;
      std::complex<double> Compile_2248 = 1./Compile_2247;
      std::complex<double> Compile_2249 = mu3US*Compile_2248;
      std::complex<double> Compile_2250 = log(Compile_2249);
      std::complex<double> Compile_2251 = 0.5 + Compile_2250;
      std::complex<double> Compile_2268 = Compile_1169 + Compile_826;
      std::complex<double> Compile_2269 = pow(Compile_2268,2);
      std::complex<double> Compile_2253 = Compile_118 + Compile_826;
      std::complex<double> Compile_2254 = \
      -0.5*(Compile_11*Compile_2253*Compile_5);
      std::complex<double> Compile_2256 = Compile_1051 + Compile_119 + \
      Compile_827;
      std::complex<double> Compile_2257 = 1./Compile_2256;
      std::complex<double> Compile_2258 = mu3US*Compile_2257;
      std::complex<double> Compile_2259 = log(Compile_2258);
      std::complex<double> Compile_2260 = 0.5 + Compile_2259;
      std::complex<double> Compile_1633 = Compile_1632*Compile_79;
      std::complex<double> Compile_1637 = -(Compile_1636*Compile_62);
      std::complex<double> Compile_1638 = Compile_1633 + Compile_1637;
      std::complex<double> Compile_1639 = 4*Compile_1638*Compile_577;
      std::complex<double> Compile_1643 = Compile_1642*Compile_62;
      std::complex<double> Compile_1647 = -(Compile_1646*Compile_79);
      std::complex<double> Compile_1648 = Compile_1643 + Compile_1647;
      std::complex<double> Compile_1659 = -4*Compile_1636*Compile_577;
      std::complex<double> Compile_1662 = 4*Compile_1632*Compile_577;
      std::complex<double> Compile_2278 = Compile_1208 + Compile_827;
      std::complex<double> Compile_2279 = 1./Compile_2278;
      std::complex<double> Compile_2280 = mu3US*Compile_2279;
      std::complex<double> Compile_2281 = log(Compile_2280);
      std::complex<double> Compile_2282 = 0.5 + Compile_2281;
      std::complex<double> Compile_811 = \
      -(Compile_140*Compile_627*Compile_80*Compile_805);
      std::complex<double> Compile_812 = -8 + Compile_811;
      std::complex<double> Compile_813 = pow(Compile_812,2);
      std::complex<double> Compile_2240 = \
      (Compile_11*Compile_1201*Compile_5*Compile_827)/16.;
      std::complex<double> Compile_2314 = Compile_167 + Compile_2190;
      std::complex<double> Compile_2315 = pow(Compile_2314,2);
      std::complex<double> Compile_857 = \
      Compile_140*Compile_627*Compile_80*Compile_805;
      std::complex<double> Compile_858 = 8 + Compile_857;
      std::complex<double> Compile_859 = pow(Compile_858,2);
      std::complex<double> Compile_2310 = Compile_1085 + Compile_167 + \
      Compile_2190;
      std::complex<double> Compile_2308 = Compile_1085 + Compile_1261 + \
      Compile_826;
      std::complex<double> Compile_2316 = Compile_168 + Compile_827;
      std::complex<double> Compile_2317 = 1./Compile_2316;
      std::complex<double> Compile_2318 = mu3US*Compile_2317;
      std::complex<double> Compile_2319 = log(Compile_2318);
      std::complex<double> Compile_2320 = 0.5 + Compile_2319;
      std::complex<double> Compile_2337 = Compile_1261 + Compile_826;
      std::complex<double> Compile_2338 = pow(Compile_2337,2);
      std::complex<double> Compile_2322 = Compile_167 + Compile_826;
      std::complex<double> Compile_2323 = \
      -0.5*(Compile_11*Compile_2322*Compile_5);
      std::complex<double> Compile_2325 = Compile_1051 + Compile_168 + \
      Compile_827;
      std::complex<double> Compile_2326 = 1./Compile_2325;
      std::complex<double> Compile_2327 = mu3US*Compile_2326;
      std::complex<double> Compile_2328 = log(Compile_2327);
      std::complex<double> Compile_2329 = 0.5 + Compile_2328;
      std::complex<double> Compile_2031 = \
      -16*m12R3dUS*lam1H3dUS*lam4H3dUS*Compile_1136;
      std::complex<double> Compile_2041 = 4*lam5H3dUS*phi2*Compile_178;
      std::complex<double> Compile_2042 = -8*m13dUS*m23dUS*lam5H3dUS*phi2;
      std::complex<double> Compile_2043 = 4*lam5H3dUS*phi2*Compile_172;
      std::complex<double> Compile_2044 = \
      8*m13dUS*lam1H3dUS*lam5H3dUS*phi2*Compile_9;
      std::complex<double> Compile_2045 = \
      -8*m23dUS*lam1H3dUS*lam5H3dUS*phi2*Compile_9;
      std::complex<double> Compile_2046 = \
      8*m13dUS*lam3H3dUS*lam5H3dUS*phi2*Compile_9;
      std::complex<double> Compile_2047 = \
      -8*m23dUS*lam3H3dUS*lam5H3dUS*phi2*Compile_9;
      std::complex<double> Compile_2052 = \
      -12*lam5H3dUS*phi2*Compile_190*Compile_191;
      std::complex<double> Compile_2053 = \
      16*lam1H3dUS*lam3H3dUS*lam5H3dUS*phi2*Compile_191;
      std::complex<double> Compile_2054 = \
      -5*lam5H3dUS*phi2*Compile_184*Compile_191;
      std::complex<double> Compile_2056 = \
      -10*lam3H3dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_191;
      std::complex<double> Compile_2073 = \
      -16*m13dUS*lam2H3dUS*lam5H3dUS*Compile_1397;
      std::complex<double> Compile_2074 = \
      16*m23dUS*lam2H3dUS*lam5H3dUS*Compile_1397;
      std::complex<double> Compile_2075 = \
      4*m13dUS*lam3H3dUS*lam5H3dUS*Compile_1397;
      std::complex<double> Compile_2076 = \
      -4*m23dUS*lam3H3dUS*lam5H3dUS*Compile_1397;
      std::complex<double> Compile_2077 = \
      4*m13dUS*lam4H3dUS*lam5H3dUS*Compile_1397;
      std::complex<double> Compile_2078 = \
      -4*m23dUS*lam4H3dUS*lam5H3dUS*Compile_1397;
      std::complex<double> Compile_2079 = \
      4*lam1H3dUS*lam3H3dUS*lam5H3dUS*Compile_1397*Compile_9;
      std::complex<double> Compile_2080 = \
      4*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_1397*Compile_9;
      std::complex<double> Compile_2081 = \
      -4*lam5H3dUS*Compile_1397*Compile_184*Compile_9;
      std::complex<double> Compile_2083 = \
      -8*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_1397*Compile_9;
      std::complex<double> Compile_2088 = \
      12*lam5H3dUS*Compile_1810*Compile_238;
      std::complex<double> Compile_2089 = \
      -8*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_1810;
      std::complex<double> Compile_2090 = lam5H3dUS*Compile_1810*Compile_184;
      std::complex<double> Compile_2092 = \
      2*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_1810;
      std::complex<double> Compile_1823 = \
      8*m12R3dUS*lam1H3dUS*phi1*Compile_608;
      std::complex<double> Compile_1824 = \
      -4*m12R3dUS*lam3H3dUS*phi1*Compile_608;
      std::complex<double> Compile_1825 = \
      2*m13dUS*lam4H3dUS*phi2*Compile_608;
      std::complex<double> Compile_1826 = \
      -2*m23dUS*lam4H3dUS*phi2*Compile_608;
      std::complex<double> Compile_1827 = \
      2*m13dUS*lam5H3dUS*phi2*Compile_608;
      std::complex<double> Compile_1828 = \
      -2*m23dUS*lam5H3dUS*phi2*Compile_608;
      std::complex<double> Compile_1829 = \
      -2*lam1H3dUS*lam4H3dUS*phi2*Compile_608*Compile_9;
      std::complex<double> Compile_1830 = \
      lam3H3dUS*lam4H3dUS*phi2*Compile_608*Compile_9;
      std::complex<double> Compile_1831 = \
      -2*lam1H3dUS*lam5H3dUS*phi2*Compile_608*Compile_9;
      std::complex<double> Compile_1832 = \
      lam3H3dUS*lam5H3dUS*phi2*Compile_608*Compile_9;
      std::complex<double> Compile_1833 = \
      -2*lam2H3dUS*lam4H3dUS*Compile_1397*Compile_608;
      std::complex<double> Compile_1834 = \
      lam3H3dUS*lam4H3dUS*Compile_1397*Compile_608;
      std::complex<double> Compile_1835 = \
      -2*lam2H3dUS*lam5H3dUS*Compile_1397*Compile_608;
      std::complex<double> Compile_1836 = \
      lam3H3dUS*lam5H3dUS*Compile_1397*Compile_608;
      std::complex<double> Compile_1847 = Compile_1846*Compile_79;
      std::complex<double> Compile_1850 = -(Compile_154*Compile_1849);
      std::complex<double> Compile_1851 = Compile_1847 + Compile_1850;
      std::complex<double> Compile_1852 = 4*Compile_1851*Compile_577;
      std::complex<double> Compile_1855 = Compile_154*Compile_1854;
      std::complex<double> Compile_1858 = -(Compile_1857*Compile_79);
      std::complex<double> Compile_1859 = Compile_1855 + Compile_1858;
      std::complex<double> Compile_1870 = -4*Compile_1849*Compile_577;
      std::complex<double> Compile_1873 = 4*Compile_1846*Compile_577;
      std::complex<double> Compile_2429 = Compile_1298 + Compile_827;
      std::complex<double> Compile_2430 = 1./Compile_2429;
      std::complex<double> Compile_2431 = mu3US*Compile_2430;
      std::complex<double> Compile_2432 = log(Compile_2431);
      std::complex<double> Compile_2433 = 0.5 + Compile_2432;
      std::complex<double> Compile_2293 = -(lam1H3dUS*phi1*Compile_788);
      std::complex<double> Compile_2294 = Compile_1670 + Compile_2293;
      std::complex<double> Compile_2296 = \
      -(Compile_21*Compile_2169*Compile_79);
      std::complex<double> Compile_2299 = -(lam3H3dUS*phi1*Compile_788);
      std::complex<double> Compile_2300 = Compile_1677 + Compile_2299;
      std::complex<double> Compile_2301 = 4*Compile_2300*Compile_79;
      std::complex<double> Compile_2460 = Compile_2190 + Compile_362;
      std::complex<double> Compile_2461 = pow(Compile_2460,2);
      std::complex<double> Compile_2456 = Compile_1040 + Compile_2190 + \
      Compile_362;
      std::complex<double> Compile_1916 = \
      (Compile_11*Compile_1201*Compile_363*Compile_6)/16.;
      std::complex<double> Compile_2453 = Compile_1040 + Compile_1311 + \
      Compile_826;
      std::complex<double> Compile_2462 = Compile_363 + Compile_827;
      std::complex<double> Compile_2463 = 1./Compile_2462;
      std::complex<double> Compile_2464 = mu3US*Compile_2463;
      std::complex<double> Compile_2465 = log(Compile_2464);
      std::complex<double> Compile_2466 = 0.5 + Compile_2465;
      std::complex<double> Compile_2483 = Compile_1311 + Compile_826;
      std::complex<double> Compile_2484 = pow(Compile_2483,2);
      std::complex<double> Compile_2468 = Compile_362 + Compile_826;
      std::complex<double> Compile_2469 = \
      -0.5*(Compile_11*Compile_2468*Compile_6);
      std::complex<double> Compile_2471 = Compile_1052 + Compile_363 + \
      Compile_827;
      std::complex<double> Compile_2472 = 1./Compile_2471;
      std::complex<double> Compile_2473 = mu3US*Compile_2472;
      std::complex<double> Compile_2474 = log(Compile_2473);
      std::complex<double> Compile_2475 = 0.5 + Compile_2474;
      std::complex<double> Compile_1930 = 4*Compile_1929*Compile_291;
      std::complex<double> Compile_1934 = -(Compile_1933*Compile_369);
      std::complex<double> Compile_1935 = Compile_1930 + Compile_1934;
      std::complex<double> Compile_1936 = 2*Compile_1935*Compile_577;
      std::complex<double> Compile_1940 = Compile_1939*Compile_369;
      std::complex<double> Compile_1944 = -2*Compile_1943*Compile_291;
      std::complex<double> Compile_1945 = Compile_1940 + Compile_1944;
      std::complex<double> Compile_1958 = -2*Compile_1933*Compile_577;
      std::complex<double> Compile_1961 = 4*Compile_1929*Compile_577;
      std::complex<double> Compile_2493 = Compile_1950 + Compile_827;
      std::complex<double> Compile_2494 = 1./Compile_2493;
      std::complex<double> Compile_2495 = mu3US*Compile_2494;
      std::complex<double> Compile_2496 = log(Compile_2495);
      std::complex<double> Compile_2497 = 0.5 + Compile_2496;
      std::complex<double> Compile_2455 = \
      (Compile_11*Compile_1201*Compile_6*Compile_827)/16.;
      std::complex<double> Compile_2528 = Compile_2190 + Compile_468;
      std::complex<double> Compile_2529 = pow(Compile_2528,2);
      std::complex<double> Compile_2524 = Compile_1040 + Compile_2190 + \
      Compile_468;
      std::complex<double> Compile_2009 = \
      (Compile_11*Compile_1201*Compile_469*Compile_6)/16.;
      std::complex<double> Compile_2522 = Compile_1040 + Compile_1414 + \
      Compile_826;
      std::complex<double> Compile_2530 = Compile_469 + Compile_827;
      std::complex<double> Compile_2531 = 1./Compile_2530;
      std::complex<double> Compile_2532 = mu3US*Compile_2531;
      std::complex<double> Compile_2533 = log(Compile_2532);
      std::complex<double> Compile_2534 = 0.5 + Compile_2533;
      std::complex<double> Compile_2551 = Compile_1414 + Compile_826;
      std::complex<double> Compile_2552 = pow(Compile_2551,2);
      std::complex<double> Compile_2536 = Compile_468 + Compile_826;
      std::complex<double> Compile_2537 = \
      -0.5*(Compile_11*Compile_2536*Compile_6);
      std::complex<double> Compile_2539 = Compile_1052 + Compile_469 + \
      Compile_827;
      std::complex<double> Compile_2540 = 1./Compile_2539;
      std::complex<double> Compile_2541 = mu3US*Compile_2540;
      std::complex<double> Compile_2542 = log(Compile_2541);
      std::complex<double> Compile_2543 = 0.5 + Compile_2542;
      std::complex<double> Compile_2020 = 8*m12R3dUS*m13dUS*lam1H3dUS*phi1;
      std::complex<double> Compile_2021 = \
      -8*m12R3dUS*m23dUS*lam1H3dUS*phi1;
      std::complex<double> Compile_2022 = \
      -4*m12R3dUS*m13dUS*lam3H3dUS*phi1;
      std::complex<double> Compile_2023 = 4*m12R3dUS*m23dUS*lam3H3dUS*phi1;
      std::complex<double> Compile_2024 = \
      -4*m12R3dUS*m13dUS*lam4H3dUS*phi1;
      std::complex<double> Compile_2025 = 4*m12R3dUS*m23dUS*lam4H3dUS*phi1;
      std::complex<double> Compile_2026 = \
      -4*m12R3dUS*m13dUS*lam5H3dUS*phi1;
      std::complex<double> Compile_2027 = 4*m12R3dUS*m23dUS*lam5H3dUS*phi1;
      std::complex<double> Compile_2028 = \
      24*m12R3dUS*Compile_1136*Compile_190;
      std::complex<double> Compile_2029 = \
      -16*m12R3dUS*lam1H3dUS*lam3H3dUS*Compile_1136;
      std::complex<double> Compile_2030 = \
      2*m12R3dUS*Compile_1136*Compile_184;
      std::complex<double> Compile_2032 = \
      4*m12R3dUS*lam3H3dUS*lam4H3dUS*Compile_1136;
      std::complex<double> Compile_2033 = \
      2*m12R3dUS*Compile_1136*Compile_217;
      std::complex<double> Compile_2034 = \
      4*m12R3dUS*lam3H3dUS*lam5H3dUS*Compile_1136;
      std::complex<double> Compile_2035 = \
      4*m12R3dUS*lam4H3dUS*lam5H3dUS*Compile_1136;
      std::complex<double> Compile_2036 = \
      -6*m12R3dUS*Compile_1136*Compile_231;
      std::complex<double> Compile_2037 = 16*lam2H3dUS*phi2*Compile_34;
      std::complex<double> Compile_2038 = -8*lam3H3dUS*phi2*Compile_34;
      std::complex<double> Compile_2039 = -8*lam4H3dUS*phi2*Compile_34;
      std::complex<double> Compile_2040 = 8*lam5H3dUS*phi2*Compile_34;
      std::complex<double> Compile_2048 = \
      8*m13dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_9;
      std::complex<double> Compile_2049 = \
      -8*m23dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_9;
      std::complex<double> Compile_2050 = \
      4*m13dUS*phi2*Compile_231*Compile_9;
      std::complex<double> Compile_2051 = \
      -4*m23dUS*phi2*Compile_231*Compile_9;
      std::complex<double> Compile_2055 = \
      16*lam1H3dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_191;
      std::complex<double> Compile_2057 = \
      -5*lam5H3dUS*phi2*Compile_191*Compile_217;
      std::complex<double> Compile_2058 = \
      4*lam1H3dUS*phi2*Compile_191*Compile_231;
      std::complex<double> Compile_2059 = 5*phi2*Compile_191*Compile_236;
      std::complex<double> Compile_2060 = \
      -24*m12R3dUS*lam1H3dUS*lam2H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2061 = \
      4*m12R3dUS*lam1H3dUS*lam3H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2062 = \
      -4*m12R3dUS*lam2H3dUS*lam3H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2063 = \
      6*m12R3dUS*phi1*Compile_10*Compile_184;
      std::complex<double> Compile_2064 = \
      4*m12R3dUS*lam1H3dUS*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2065 = \
      -4*m12R3dUS*lam2H3dUS*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2066 = \
      12*m12R3dUS*lam3H3dUS*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2067 = \
      6*m12R3dUS*phi1*Compile_10*Compile_217;
      std::complex<double> Compile_2068 = \
      4*m12R3dUS*lam1H3dUS*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2069 = \
      -36*m12R3dUS*lam2H3dUS*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2070 = \
      4*m12R3dUS*lam3H3dUS*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2071 = \
      4*m12R3dUS*lam4H3dUS*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2072 = \
      -10*m12R3dUS*phi1*Compile_10*Compile_231;
      std::complex<double> Compile_2082 = \
      4*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_1397*Compile_9;
      std::complex<double> Compile_2084 = \
      -4*lam5H3dUS*Compile_1397*Compile_217*Compile_9;
      std::complex<double> Compile_2085 = \
      -8*lam1H3dUS*Compile_1397*Compile_231*Compile_9;
      std::complex<double> Compile_2086 = \
      16*lam2H3dUS*Compile_1397*Compile_231*Compile_9;
      std::complex<double> Compile_2087 = \
      4*Compile_1397*Compile_236*Compile_9;
      std::complex<double> Compile_2091 = \
      -8*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_1810;
      std::complex<double> Compile_2093 = lam5H3dUS*Compile_1810*Compile_217;
      std::complex<double> Compile_2094 = \
      4*lam2H3dUS*Compile_1810*Compile_231;
      std::complex<double> Compile_2095 = -(Compile_1810*Compile_236);
      std::complex<double> Compile_2120 = 4*Compile_2119*Compile_291;
      std::complex<double> Compile_2123 = -(Compile_2122*Compile_431);
      std::complex<double> Compile_2124 = Compile_2120 + Compile_2123;
      std::complex<double> Compile_2125 = 2*Compile_2124*Compile_577;
      std::complex<double> Compile_2128 = Compile_2127*Compile_431;
      std::complex<double> Compile_2131 = -2*Compile_2130*Compile_291;
      std::complex<double> Compile_2132 = Compile_2128 + Compile_2131;
      std::complex<double> Compile_2145 = -2*Compile_2122*Compile_577;
      std::complex<double> Compile_2148 = 4*Compile_2119*Compile_577;
      std::complex<double> Compile_2577 = Compile_2137 + Compile_827;
      std::complex<double> Compile_2578 = 1./Compile_2577;
      std::complex<double> Compile_2579 = mu3US*Compile_2578;
      std::complex<double> Compile_2580 = log(Compile_2579);
      std::complex<double> Compile_2581 = 0.5 + Compile_2580;
      std::complex<double> Compile_2508 = Compile_1969 + Compile_2293;
      std::complex<double> Compile_2510 = \
      -2*lam5H3dUS*Compile_2169*Compile_291;
      std::complex<double> Compile_2513 = -(phi1*Compile_180*Compile_788);
      std::complex<double> Compile_2514 = Compile_1677 + Compile_2513;
      std::complex<double> Compile_2515 = 4*Compile_2514*Compile_291;
      std::complex<double> Compile_982 = lam4H3dUS*Compile_10;
      std::complex<double> Compile_1507 = \
      -2*Compile_1506*Compile_577*Compile_71;
      std::complex<double> Compile_1509 = -3*lam1H3dUS*phi1*Compile_661;
      std::complex<double> Compile_1510 = Compile_1508 + Compile_1509;
      std::complex<double> Compile_1511 = Compile_1510*Compile_661;
      std::complex<double> Compile_1512 = Compile_1507 + Compile_1511;
      std::complex<double> Compile_1514 = \
      -(Compile_1506*Compile_661*Compile_71);
      std::complex<double> Compile_1516 = -(phi1*Compile_661*Compile_71);
      std::complex<double> Compile_1517 = Compile_1515 + Compile_1516;
      std::complex<double> Compile_1518 = 4*Compile_1517*Compile_577;
      std::complex<double> Compile_1519 = Compile_1514 + Compile_1518;
      std::complex<double> Compile_1520 = 2*Compile_1519*Compile_577;
      std::complex<double> Compile_2643 = 2*Compile_655;
      std::complex<double> Compile_2644 = Compile_2643 + Compile_827;
      std::complex<double> Compile_2645 = 1./Compile_2644;
      std::complex<double> Compile_2646 = mu3US*Compile_2645;
      std::complex<double> Compile_2647 = log(Compile_2646);
      std::complex<double> Compile_2648 = 0.5 + Compile_2647;
      std::complex<double> Compile_2601 = phi1*phi2*Compile_71;
      std::complex<double> Compile_2602 = Compile_2601 + Compile_297;
      std::complex<double> Compile_2603 = \
      16*m13dUS*phi1*Compile_2602*Compile_71;
      std::complex<double> Compile_2604 = -2*m23dUS*phi1*Compile_71;
      std::complex<double> Compile_2605 = -(Compile_1136*Compile_184);
      std::complex<double> Compile_2606 = \
      -2*lam3H3dUS*lam4H3dUS*Compile_1136;
      std::complex<double> Compile_2607 = -(Compile_1136*Compile_217);
      std::complex<double> Compile_2608 = \
      -2*lam3H3dUS*lam5H3dUS*Compile_1136;
      std::complex<double> Compile_2609 = \
      -2*lam4H3dUS*lam5H3dUS*Compile_1136;
      std::complex<double> Compile_2610 = -(Compile_1136*Compile_231);
      std::complex<double> Compile_2611 = \
      6*lam1H3dUS*Compile_1136*Compile_71;
      std::complex<double> Compile_2612 = -12*m12R3dUS*lam2H3dUS*phi2;
      std::complex<double> Compile_2613 = 2*m12R3dUS*lam3H3dUS*phi2;
      std::complex<double> Compile_2614 = 2*m12R3dUS*lam4H3dUS*phi2;
      std::complex<double> Compile_2615 = 2*m12R3dUS*lam5H3dUS*phi2;
      std::complex<double> Compile_2616 = \
      6*lam2H3dUS*lam3H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2617 = -(phi1*Compile_10*Compile_184);
      std::complex<double> Compile_2618 = \
      6*lam2H3dUS*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2619 = \
      -2*lam3H3dUS*lam4H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2620 = -(phi1*Compile_10*Compile_217);
      std::complex<double> Compile_2621 = \
      6*lam2H3dUS*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2622 = \
      -2*lam3H3dUS*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2623 = \
      -2*lam4H3dUS*lam5H3dUS*phi1*Compile_10;
      std::complex<double> Compile_2624 = -(phi1*Compile_10*Compile_231);
      std::complex<double> Compile_2625 = Compile_2604 + Compile_2605 + \
      Compile_2606 + Compile_2607 + Compile_2608 + Compile_2609 + \
      Compile_2610 + Compile_2611 + Compile_2612 + Compile_2613 + \
      Compile_2614 + Compile_2615 + Compile_2616 + Compile_2617 + \
      Compile_2618 + Compile_2619 + Compile_2620 + Compile_2621 + \
      Compile_2622 + Compile_2623 + Compile_2624;
      std::complex<double> Compile_2626 = -8*Compile_2625*Compile_577;
      std::complex<double> Compile_2627 = Compile_2603 + Compile_2626;
      std::complex<double> Compile_2628 = 2*Compile_2627*Compile_577;
      std::complex<double> Compile_2629 = 8*phi1*Compile_34*Compile_587;
      std::complex<double> Compile_2630 = pow(Compile_71,2);
      std::complex<double> Compile_2631 = Compile_28 + Compile_29 + \
      Compile_33 + Compile_568 + Compile_571 + Compile_583 + Compile_584 + \
      Compile_658 + Compile_97 + Compile_982;
      std::complex<double> Compile_2632 = \
      -4*phi1*Compile_10*Compile_2630*Compile_2631;
      std::complex<double> Compile_2633 = -18*lam1H3dUS*Compile_9;
      std::complex<double> Compile_2634 = 3*lam3H3dUS*Compile_9;
      std::complex<double> Compile_2635 = 3*lam4H3dUS*Compile_9;
      std::complex<double> Compile_2636 = 3*lam5H3dUS*Compile_9;
      std::complex<double> Compile_2637 = Compile_2633 + Compile_2634 + \
      Compile_2635 + Compile_2636 + Compile_28 + Compile_29 + Compile_33 + \
      Compile_571 + Compile_658 + Compile_982;
      std::complex<double> Compile_2638 = \
      4*m12R3dUS*phi2*Compile_2637*Compile_71;
      std::complex<double> Compile_2639 = Compile_2629 + Compile_2632 + \
      Compile_2638;
      std::complex<double> Compile_2164 = -3*lam1H3dUS*phi1*Compile_788;
      std::complex<double> Compile_2165 = Compile_1508 + Compile_2164;
      std::complex<double> Compile_2166 = Compile_2165*Compile_788;
      std::complex<double> Compile_2170 = \
      -2*Compile_2169*Compile_577*Compile_71;
      std::complex<double> Compile_2171 = Compile_2166 + Compile_2170;
      std::complex<double> Compile_2173 = -(phi1*Compile_71*Compile_788);
      std::complex<double> Compile_2174 = Compile_1515 + Compile_2173;
      std::complex<double> Compile_2175 = 4*Compile_2174*Compile_577;
      std::complex<double> Compile_2176 = \
      -(Compile_2169*Compile_71*Compile_788);
      std::complex<double> Compile_2177 = Compile_2175 + Compile_2176;
      std::complex<double> Compile_2178 = 2*Compile_2177*Compile_577;
      std::complex<double> Compile_2657 = 2*Compile_827;
      std::complex<double> Compile_2658 = Compile_2657 + Compile_655;
      std::complex<double> Compile_2659 = 1./Compile_2658;
      std::complex<double> Compile_2660 = mu3US*Compile_2659;
      std::complex<double> Compile_2661 = log(Compile_2660);
      std::complex<double> Compile_2662 = 0.5 + Compile_2661;
      
      return (8*Compile_1011*Compile_1012*(Compile_1018 + (Compile_1024 + \
      Compile_1029 + Compile_1034 + Compile_1038 + \
      Compile_1042*Compile_1059 + Compile_1059*Compile_1088 + (Compile_1042 \
      + Compile_1044 + Compile_1091)*Compile_1104)/3.))/3. + \
      (8*Compile_1011*Compile_1012*(Compile_1018 + (Compile_1024 + \
      Compile_1029 + Compile_1034 + Compile_1038 + \
      2*Compile_1059*Compile_1088 + (Compile_1044 + Compile_1088 + \
      Compile_1091)*Compile_1104)/3.))/3. + \
      (Compile_1011*Compile_113*Compile_1133*Compile_116*Compile_120*\
      Compile_1215)/4. + \
      (Compile_1011*Compile_1133*Compile_116*Compile_120*Compile_1215*\
      Compile_124)/4. + \
      (Compile_1011*Compile_113*Compile_1133*Compile_1305*Compile_165*\
      Compile_169)/4. + \
      (Compile_1011*Compile_1133*Compile_124*Compile_1305*Compile_165*\
      Compile_169)/4. + \
      (Compile_1133*Compile_117*Compile_1333*Compile_347*Compile_348)/4. + \
      (Compile_1133*Compile_117*Compile_1346*Compile_347*Compile_348)/4. + \
      (Compile_1133*Compile_117*Compile_1333*Compile_348*Compile_352)/4. + \
      (Compile_1133*Compile_117*Compile_1346*Compile_348*Compile_352)/4. + \
      (Compile_1133*Compile_1372*Compile_166*Compile_348*Compile_356)/4. + \
      (Compile_1133*Compile_1385*Compile_166*Compile_348*Compile_356)/4. + \
      (Compile_1133*Compile_1372*Compile_166*Compile_348*Compile_360)/4. + \
      (Compile_1133*Compile_1385*Compile_166*Compile_348*Compile_360)/4. + \
      (Compile_1133*Compile_117*Compile_1436*Compile_453*Compile_454)/4. + \
      (Compile_1133*Compile_117*Compile_1448*Compile_453*Compile_454)/4. + \
      (Compile_1133*Compile_117*Compile_1436*Compile_454*Compile_458)/4. + \
      (Compile_1133*Compile_117*Compile_1448*Compile_454*Compile_458)/4. + \
      (Compile_1133*Compile_1474*Compile_166*Compile_454*Compile_462)/4. + \
      (Compile_1133*Compile_1486*Compile_166*Compile_454*Compile_462)/4. + \
      (Compile_1133*Compile_1474*Compile_166*Compile_454*Compile_466)/4. + \
      (Compile_1133*Compile_1486*Compile_166*Compile_454*Compile_466)/4. + \
      (Compile_1129*Compile_116*Compile_120*Compile_2*Compile_5*Compile_7)/\
      4. + (Compile_1229*Compile_165*Compile_169*Compile_2*Compile_5*\
      Compile_7)/4. + (Compile_1229*Compile_169*Compile_2*Compile_5*pow(-4 \
      - Compile_141*Compile_64,2)*Compile_7)/4. + \
      (Compile_1129*Compile_120*Compile_2*Compile_5*pow(-4 - \
      Compile_64*Compile_66,2)*Compile_7)/4. + \
      (Compile_1009*Compile_11*Compile_3*Compile_8)/256. + \
      (Compile_117*Compile_119*Compile_16*Compile_347*Compile_348*Compile_5*\
      Compile_8)/128. + \
      (Compile_117*Compile_119*Compile_16*Compile_348*Compile_352*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_168*Compile_348*Compile_356*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_168*Compile_348*Compile_360*Compile_5*\
      Compile_8)/128. + \
      (Compile_117*Compile_16*Compile_347*Compile_348*Compile_363*Compile_5*\
      Compile_8)/128. + \
      (Compile_117*Compile_16*Compile_348*Compile_352*Compile_363*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_348*Compile_356*Compile_363*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_348*Compile_360*Compile_363*Compile_5*\
      Compile_8)/128. + \
      (Compile_117*Compile_119*Compile_16*Compile_453*Compile_454*Compile_5*\
      Compile_8)/128. + \
      (Compile_117*Compile_119*Compile_16*Compile_454*Compile_458*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_168*Compile_454*Compile_462*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_168*Compile_454*Compile_466*Compile_5*\
      Compile_8)/128. + \
      (Compile_117*Compile_16*Compile_453*Compile_454*Compile_469*Compile_5*\
      Compile_8)/128. + \
      (Compile_117*Compile_16*Compile_454*Compile_458*Compile_469*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_454*Compile_462*Compile_469*Compile_5*\
      Compile_8)/128. + \
      (Compile_16*Compile_166*Compile_454*Compile_466*Compile_469*Compile_5*\
      Compile_8)/128. - \
      (Compile_1949*Compile_1955*Compile_292*pow(2*Compile_291*( Compile_1977 - lam5H3dUS*Compile_1506*Compile_369) + \
      Compile_369*(Compile_1968 + \
      Compile_1970*Compile_369),2)*Compile_611*Compile_637*Compile_8)/48. - \
      (Compile_2136*Compile_2142*Compile_292*pow(2*Compile_291*( Compile_1977 - lam5H3dUS*Compile_1506*Compile_431) + \
      Compile_431*(Compile_1968 + \
      Compile_1970*Compile_431),2)*Compile_611*Compile_637*Compile_8)/48. + \
      (Compile_117*Compile_119*Compile_16*Compile_5*Compile_636*Compile_637*\
      Compile_8)/64. + \
      (Compile_16*Compile_166*Compile_168*Compile_5*Compile_637*Compile_641*\
      Compile_8)/64. + \
      (Compile_18*Compile_348*Compile_363*Compile_6*Compile_637*Compile_645*\
      Compile_8)/128. + \
      (Compile_18*Compile_454*Compile_469*Compile_6*Compile_637*Compile_649*\
      Compile_8)/128. + \
      (Compile_18*Compile_348*Compile_6*Compile_637*Compile_645*Compile_655*\
      Compile_8)/128. + \
      (Compile_18*Compile_454*Compile_6*Compile_637*Compile_649*Compile_655*\
      Compile_8)/128. + \
      (Compile_117*Compile_16*Compile_5*Compile_637*Compile_653*Compile_655*\
      Compile_8)/64. - \
      (Compile_120*Compile_1656*Compile_24*Compile_611*Compile_637*pow(\
      Compile_1639 + Compile_1648*Compile_661,2)*Compile_8)/96. - \
      (Compile_169*Compile_1867*Compile_24*Compile_611*Compile_637*pow(\
      Compile_1852 + Compile_1859*Compile_661,2)*Compile_8)/96. - \
      (Compile_1949*Compile_1955*Compile_292*Compile_611*Compile_637*pow(\
      Compile_1936 + Compile_1945*Compile_661,2)*Compile_8)/48. - \
      (Compile_2136*Compile_2142*Compile_292*Compile_611*Compile_637*pow(\
      Compile_2125 + Compile_2132*Compile_661,2)*Compile_8)/48. - \
      (Compile_1498*Compile_2194*Compile_2662*Compile_637*pow(Compile_2178 \
      + Compile_2171*Compile_661,2)*Compile_8)/48. - \
      (Compile_1949*Compile_1955*Compile_292*Compile_611*Compile_637*pow(\
      Compile_369*(Compile_1958 + Compile_1939*Compile_661) + \
      2*Compile_291*(Compile_1961 - \
      Compile_1943*Compile_661),2)*Compile_8)/48. - \
      (Compile_2136*Compile_2142*Compile_292*Compile_611*Compile_637*pow(\
      Compile_431*(Compile_2145 + Compile_2127*Compile_661) + \
      2*Compile_291*(Compile_2148 - \
      Compile_2130*Compile_661),2)*Compile_8)/48. + \
      (Compile_16*Compile_166*Compile_5*Compile_637*Compile_655*Compile_696*\
      Compile_8)/64. + (Compile_11*Compile_13*Compile_7*Compile_8)/48. + \
      (Compile_113*Compile_116*Compile_119*Compile_120*Compile_18*Compile_7*\
      Compile_8)/128. + \
      (Compile_116*Compile_119*Compile_120*Compile_124*Compile_18*Compile_7*\
      Compile_8)/128. + \
      (Compile_113*Compile_165*Compile_168*Compile_169*Compile_18*Compile_7*\
      Compile_8)/128. + \
      (Compile_124*Compile_165*Compile_168*Compile_169*Compile_18*Compile_7*\
      Compile_8)/128. + \
      (Compile_16*Compile_18*Compile_3*Compile_7*Compile_8)/24. + \
      (Compile_11*Compile_2*Compile_3*Compile_7*Compile_8)/48. + \
      (Compile_292*Compile_330*(3*lam1H3dUS*pow(Compile_323,4) + \
      Compile_332 + \
      48*Compile_328*Compile_333*Compile_71)*Compile_8)/(64.*pow( Compile_341,3)) + (Compile_292*Compile_435*(Compile_332 + \
      3*lam1H3dUS*pow(Compile_431,4) + \
      48*Compile_333*Compile_433*Compile_71)*Compile_8)/(64.*pow( Compile_447,3)) + \
      (Compile_578*Compile_619*Compile_632*(3*lam1H3dUS*pow(Compile_616,4) \
      + Compile_621 + 48*Compile_617*Compile_622*Compile_71)*Compile_8)/64. \
      + (Compile_143*Compile_148*Compile_24*(3*lam1H3dUS*pow(Compile_126,4) \
      + Compile_70 + 12*Compile_127*Compile_71*Compile_72)*Compile_8)/32. + \
      (Compile_325*Compile_348*Compile_363*Compile_611*Compile_637*Compile_655*(Compile_661*(8*Compile_291*Compile_729 + \
      Compile_369*Compile_733) + 8*Compile_577*(Compile_369*Compile_739 + \
      Compile_291*Compile_743))*Compile_8)/64. + \
      (Compile_325*Compile_348*Compile_363*Compile_611*Compile_637*Compile_655*(Compile_369*(Compile_661*Compile_733 + \
      8*Compile_577*Compile_739) + 8*Compile_291*(Compile_661*Compile_729 + \
      Compile_577*Compile_743))*Compile_8)/64. + \
      (Compile_325*Compile_454*Compile_469*Compile_611*Compile_637*Compile_655*(Compile_661*(8*Compile_291*Compile_761 + \
      Compile_431*Compile_764) + 8*Compile_577*(Compile_431*Compile_769 + \
      Compile_291*Compile_772))*Compile_8)/64. + \
      (Compile_325*Compile_454*Compile_469*Compile_611*Compile_637*Compile_655*(Compile_431*(Compile_661*Compile_764 + \
      8*Compile_577*Compile_769) + 8*Compile_291*(Compile_661*Compile_761 + \
      Compile_577*Compile_772))*Compile_8)/64. - \
      (Compile_1498*Compile_2194*Compile_2662*Compile_637*pow(Compile_2628 \
      + Compile_2639*Compile_788,2)*Compile_8)/24. - \
      (Compile_169*Compile_1867*Compile_24*Compile_611*Compile_637*pow(\
      Compile_154*(Compile_1669 + 2*Compile_154*Compile_1672) + \
      (Compile_1680 - \
      Compile_1506*Compile_154*Compile_21)*Compile_79,2)*Compile_8)/96. + \
      (Compile_117*Compile_119*Compile_325*Compile_348*Compile_363*Compile_64*((2*Compile_291*Compile_371 + Compile_369*Compile_375)*Compile_62 \
      + 2*(Compile_369*Compile_381 + \
      Compile_291*Compile_385)*Compile_79)*Compile_8)/32. + \
      (Compile_166*Compile_168*Compile_325*Compile_348*Compile_363*Compile_64*(Compile_154*(2*Compile_291*Compile_404 + Compile_369*Compile_407) \
      + 2*(Compile_369*Compile_412 + \
      Compile_291*Compile_415)*Compile_79)*Compile_8)/32. + \
      (Compile_117*Compile_119*Compile_325*Compile_454*Compile_469*Compile_64*((2*Compile_291*Compile_473 + Compile_431*Compile_476)*Compile_62 \
      + 2*(Compile_431*Compile_481 + \
      Compile_291*Compile_484)*Compile_79)*Compile_8)/32. + \
      (Compile_166*Compile_168*Compile_325*Compile_454*Compile_469*Compile_64*(Compile_154*(2*Compile_291*Compile_502 + Compile_431*Compile_505) \
      + 2*(Compile_431*Compile_509 + \
      Compile_291*Compile_512)*Compile_79)*Compile_8)/32. - \
      (Compile_120*Compile_1656*Compile_24*Compile_611*Compile_637*pow(\
      Compile_62*(Compile_1669 + 2*Compile_1672*Compile_62) + (Compile_1680 \
      - Compile_1506*Compile_21*Compile_62)*Compile_79,2)*Compile_8)/96. - \
      (Compile_120*Compile_1656*Compile_24*Compile_611*Compile_637*pow(\
      Compile_62*(Compile_1659 + Compile_1642*Compile_661) + (Compile_1662 \
      - Compile_1646*Compile_661)*Compile_79,2)*Compile_8)/96. - \
      (Compile_169*Compile_1867*Compile_24*Compile_611*Compile_637*pow(\
      Compile_154*(Compile_1870 + Compile_1854*Compile_661) + (Compile_1873 \
      - Compile_1857*Compile_661)*Compile_79,2)*Compile_8)/96. + \
      (Compile_117*Compile_119*Compile_611*Compile_637*Compile_64*Compile_655*(Compile_62*(Compile_661*Compile_668 + 2*Compile_577*Compile_674) \
      + 2*(Compile_661*Compile_664 + \
      Compile_577*Compile_678)*Compile_79)*Compile_8)/32. + \
      (Compile_166*Compile_168*Compile_611*Compile_637*Compile_64*Compile_655*(Compile_154*(Compile_661*Compile_702 + \
      2*Compile_577*Compile_707) + 2*(Compile_661*Compile_699 + \
      Compile_577*Compile_710)*Compile_79)*Compile_8)/32. + \
      (Compile_117*Compile_119*Compile_325*Compile_348*Compile_363*Compile_64*(Compile_369*(Compile_375*Compile_62 + 2*Compile_381*Compile_79) + \
      2*Compile_291*(Compile_371*Compile_62 + \
      Compile_385*Compile_79))*Compile_8)/32. + \
      (Compile_166*Compile_168*Compile_325*Compile_348*Compile_363*Compile_64*(Compile_369*(Compile_154*Compile_407 + 2*Compile_412*Compile_79) \
      + 2*Compile_291*(Compile_154*Compile_404 + \
      Compile_415*Compile_79))*Compile_8)/32. + \
      (Compile_117*Compile_119*Compile_325*Compile_454*Compile_469*Compile_64*(Compile_431*(Compile_476*Compile_62 + 2*Compile_481*Compile_79) + \
      2*Compile_291*(Compile_473*Compile_62 + \
      Compile_484*Compile_79))*Compile_8)/32. + \
      (Compile_166*Compile_168*Compile_325*Compile_454*Compile_469*Compile_64*(Compile_431*(Compile_154*Compile_505 + 2*Compile_509*Compile_79) \
      + 2*Compile_291*(Compile_154*Compile_502 + \
      Compile_512*Compile_79))*Compile_8)/32. + \
      (Compile_117*Compile_119*Compile_611*Compile_637*Compile_64*Compile_655*(Compile_661*(Compile_62*Compile_668 + 2*Compile_664*Compile_79) \
      + 2*Compile_577*(Compile_62*Compile_674 + \
      Compile_678*Compile_79))*Compile_8)/32. + \
      (Compile_166*Compile_168*Compile_611*Compile_637*Compile_64*Compile_655*(Compile_661*(Compile_154*Compile_702 + 2*Compile_699*Compile_79) \
      + 2*Compile_577*(Compile_154*Compile_707 + \
      Compile_710*Compile_79))*Compile_8)/32. + \
      (Compile_1133*Compile_117*Compile_637*Compile_653*((Compile_1201*(\
      Compile_1337 - \
      (Compile_1201*Compile_1593*Compile_16)/8.)*Compile_655)/4. - \
      (Compile_119*Compile_1596*Compile_16*Compile_8)/32. + \
      (Compile_1606*Compile_1624*Compile_8)/16. - \
      (Compile_1615*(Compile_1044 + Compile_1609 + \
      Compile_1624)*Compile_8)/16.))/2. + \
      (Compile_1133*Compile_166*Compile_637*Compile_696*((Compile_1201*(\
      Compile_1376 - \
      (Compile_1201*Compile_16*Compile_1686)/8.)*Compile_655)/4. - \
      (Compile_16*Compile_168*Compile_1688*Compile_8)/32. + \
      (Compile_1698*Compile_1716*Compile_8)/16. - \
      (Compile_1707*(Compile_1044 + Compile_1701 + \
      Compile_1716)*Compile_8)/16.))/2. + \
      (Compile_1133*Compile_348*Compile_637*Compile_645*((Compile_1201*(-0.125*(Compile_1201*Compile_18*Compile_1889) + \
      Compile_1916)*Compile_655)/4. + \
      (Compile_1902*Compile_1921*Compile_8)/16. - \
      (Compile_1911*(Compile_1046 + Compile_1905 + \
      Compile_1921)*Compile_8)/16. - \
      (Compile_18*Compile_1892*Compile_363*Compile_8)/32.))/4. + \
      (Compile_1133*Compile_454*Compile_637*Compile_649*((Compile_1201*(-0.125*(Compile_1201*Compile_18*Compile_1983) + \
      Compile_2009)*Compile_655)/4. + \
      (Compile_1995*Compile_2014*Compile_8)/16. - \
      (Compile_2004*(Compile_1046 + Compile_1998 + \
      Compile_2014)*Compile_8)/16. - \
      (Compile_18*Compile_1985*Compile_469*Compile_8)/32.))/4. + \
      (Compile_1133*Compile_117*Compile_636*Compile_637*((Compile_119*\
      Compile_1201*(Compile_1595 - \
      (Compile_1201*Compile_1596*Compile_16)/8.))/4. + \
      (Compile_1601*Compile_1606*Compile_8)/16. - ((Compile_1044 + \
      Compile_1601 + Compile_1609)*Compile_1615*Compile_8)/16. - \
      (Compile_1593*Compile_16*Compile_655*Compile_8)/32.))/2. + \
      (Compile_1133*Compile_166*Compile_637*Compile_641*((Compile_1201*\
      Compile_168*(Compile_1595 - \
      (Compile_1201*Compile_16*Compile_1688)/8.))/4. + \
      (Compile_1693*Compile_1698*Compile_8)/16. - ((Compile_1044 + \
      Compile_1693 + Compile_1701)*Compile_1707*Compile_8)/16. - \
      (Compile_16*Compile_1686*Compile_655*Compile_8)/32.))/2. + \
      (Compile_1133*Compile_348*Compile_637*Compile_645*((Compile_1201*(\
      Compile_1891 - \
      (Compile_1201*Compile_18*Compile_1892)/8.)*Compile_363)/4. + \
      (Compile_1897*Compile_1902*Compile_8)/16. - ((Compile_1046 + \
      Compile_1897 + Compile_1905)*Compile_1911*Compile_8)/16. - \
      (Compile_18*Compile_1889*Compile_655*Compile_8)/32.))/4. + \
      (Compile_1133*Compile_454*Compile_637*Compile_649*((Compile_1201*(\
      Compile_1891 - \
      (Compile_1201*Compile_18*Compile_1985)/8.)*Compile_469)/4. + \
      (Compile_1990*Compile_1995*Compile_8)/16. - ((Compile_1046 + \
      Compile_1990 + Compile_1998)*Compile_2004*Compile_8)/16. - \
      (Compile_18*Compile_1983*Compile_655*Compile_8)/32.))/4. + \
      (Compile_1132*Compile_1133*Compile_1242*Compile_166*Compile_2*Compile_5*Compile_7*(Compile_1147 + \
      2*(-0.0625*(Compile_1246*Compile_143*Compile_166*Compile_8) + \
      (Compile_1253*(Compile_1153 + Compile_167)*Compile_8)/16. + \
      (Compile_16*Compile_168*Compile_8)/32.)))/4. + \
      (Compile_1132*Compile_1133*Compile_1146*Compile_117*Compile_2*Compile_5*Compile_7*(Compile_1147 + 2*((Compile_1159*(Compile_1153 + \
      Compile_118)*Compile_8)/16. + (Compile_119*Compile_16*Compile_8)/32. \
      - (Compile_1151*Compile_117*Compile_68*Compile_8)/16.)))/4. + \
      (8*Compile_1011*Compile_1012*(Compile_1018 + (Compile_1024 + \
      Compile_1029 + Compile_1034 + Compile_1038 + \
      2*Compile_1042*Compile_1059 + (Compile_1046 - \
      (Compile_1019*Compile_5*Compile_6)/4.)*((-3*Compile_1069*Compile_1077*\
      Compile_8)/16. - ((-2*Compile_1069 + 3*(-0.5*(Compile_1019*Compile_3) \
      - Compile_1019*Compile_5*Compile_6))*Compile_8)/64.))/3.))/3. + \
      (Compile_578*(Compile_621 + 3*lam1H3dUS*pow(Compile_788,4) + \
      48*Compile_622*Compile_71*Compile_790)*Compile_792*Compile_8*Compile_803)/64. - \
      (Compile_1949*Compile_2497*Compile_292*pow(2*Compile_291*( Compile_2515 - lam5H3dUS*Compile_2169*Compile_369) + \
      Compile_369*(Compile_2510 + \
      Compile_2508*Compile_369),2)*Compile_611*Compile_8*Compile_809)/48. - \
      (Compile_2136*Compile_2581*Compile_292*pow(2*Compile_291*( Compile_2515 - lam5H3dUS*Compile_2169*Compile_431) + \
      Compile_431*(Compile_2510 + \
      Compile_2508*Compile_431),2)*Compile_611*Compile_8*Compile_809)/48. - \
      (Compile_1498*Compile_1545*Compile_2648*pow(Compile_2628 + \
      Compile_2639*Compile_661,2)*Compile_8*Compile_809)/24. - \
      (Compile_1498*Compile_1545*Compile_2648*pow(Compile_1520 + \
      Compile_1512*Compile_788,2)*Compile_8*Compile_809)/48. - \
      (Compile_120*Compile_2282*Compile_24*Compile_611*pow(Compile_1639 + \
      Compile_1648*Compile_788,2)*Compile_8*Compile_809)/96. - \
      (Compile_169*Compile_24*Compile_2433*Compile_611*pow(Compile_1852 + \
      Compile_1859*Compile_788,2)*Compile_8*Compile_809)/96. - \
      (Compile_1949*Compile_2497*Compile_292*Compile_611*pow(Compile_1936 + \
      Compile_1945*Compile_788,2)*Compile_8*Compile_809)/48. - \
      (Compile_2136*Compile_2581*Compile_292*Compile_611*pow(Compile_2125 + \
      Compile_2132*Compile_788,2)*Compile_8*Compile_809)/48. - \
      (Compile_1949*Compile_2497*Compile_292*Compile_611*pow(Compile_369*(\
      Compile_1958 + Compile_1939*Compile_788) + \
      2*Compile_291*(Compile_1961 - \
      Compile_1943*Compile_788),2)*Compile_8*Compile_809)/48. - \
      (Compile_2136*Compile_2581*Compile_292*Compile_611*pow(Compile_431*(\
      Compile_2145 + Compile_2127*Compile_788) + \
      2*Compile_291*(Compile_2148 - \
      Compile_2130*Compile_788),2)*Compile_8*Compile_809)/48. - \
      (Compile_169*Compile_24*Compile_2433*Compile_611*pow(Compile_154*(2*\
      Compile_154*Compile_2294 + Compile_2296) + \
      (-(Compile_154*Compile_21*Compile_2169) + \
      Compile_2301)*Compile_79,2)*Compile_8*Compile_809)/96. - \
      (Compile_120*Compile_2282*Compile_24*Compile_611*pow(Compile_62*(\
      Compile_2296 + 2*Compile_2294*Compile_62) + (Compile_2301 - \
      Compile_21*Compile_2169*Compile_62)*Compile_79,2)*Compile_8*Compile_809)/96. - \
      (Compile_120*Compile_2282*Compile_24*Compile_611*pow(Compile_62*(\
      Compile_1659 + Compile_1642*Compile_788) + (Compile_1662 - \
      Compile_1646*Compile_788)*Compile_79,2)*Compile_8*Compile_809)/96. - \
      (Compile_169*Compile_24*Compile_2433*Compile_611*pow(Compile_154*(\
      Compile_1870 + Compile_1854*Compile_788) + (Compile_1873 - \
      Compile_1857*Compile_788)*Compile_79,2)*Compile_8*Compile_809)/96. + \
      (Compile_117*Compile_119*Compile_16*Compile_5*Compile_8*Compile_808*\
      Compile_809)/64. + \
      (Compile_16*Compile_166*Compile_168*Compile_5*Compile_8*Compile_809*\
      Compile_813)/64. + \
      (Compile_18*Compile_348*Compile_363*Compile_6*Compile_8*Compile_809*\
      Compile_817)/128. + \
      (Compile_18*Compile_454*Compile_469*Compile_6*Compile_8*Compile_809*\
      Compile_821)/128. + \
      (Compile_18*Compile_348*Compile_6*Compile_8*Compile_809*Compile_817*\
      Compile_827)/128. + \
      (Compile_18*Compile_454*Compile_6*Compile_8*Compile_809*Compile_821*\
      Compile_827)/128. + \
      (Compile_117*Compile_16*Compile_5*Compile_8*Compile_809*Compile_825*\
      Compile_827)/64. + \
      (Compile_1133*Compile_117*Compile_809*Compile_825*(-0.03125*( Compile_119*Compile_16*Compile_2241*Compile_8) + \
      (Compile_2251*Compile_2269*Compile_8)/16. - \
      (Compile_2260*(Compile_1044 + Compile_2254 + \
      Compile_2269)*Compile_8)/16. + (Compile_1201*(Compile_1337 - \
      (Compile_1201*Compile_16*Compile_2238)/8.)*Compile_827)/4.))/2. + \
      (Compile_1133*Compile_348*Compile_809*Compile_817*((Compile_2466*\
      Compile_2484*Compile_8)/16. - (Compile_2475*(Compile_1046 + \
      Compile_2469 + Compile_2484)*Compile_8)/16. - \
      (Compile_18*Compile_2456*Compile_363*Compile_8)/32. + \
      (Compile_1201*(Compile_1916 - \
      (Compile_1201*Compile_18*Compile_2453)/8.)*Compile_827)/4.))/4. + \
      (Compile_1133*Compile_454*Compile_809*Compile_821*((Compile_2534*\
      Compile_2552*Compile_8)/16. - (Compile_2543*(Compile_1046 + \
      Compile_2537 + Compile_2552)*Compile_8)/16. - \
      (Compile_18*Compile_2524*Compile_469*Compile_8)/32. + \
      (Compile_1201*(Compile_2009 - \
      (Compile_1201*Compile_18*Compile_2522)/8.)*Compile_827)/4.))/4. + \
      (Compile_1133*Compile_117*Compile_808*Compile_809*((Compile_119*\
      Compile_1201*(Compile_2240 - \
      (Compile_1201*Compile_16*Compile_2241)/8.))/4. + \
      (Compile_2246*Compile_2251*Compile_8)/16. - ((Compile_1044 + \
      Compile_2246 + Compile_2254)*Compile_2260*Compile_8)/16. - \
      (Compile_16*Compile_2238*Compile_8*Compile_827)/32.))/2. + \
      (Compile_1133*Compile_166*Compile_809*Compile_813*((Compile_1201*\
      Compile_168*(Compile_2240 - \
      (Compile_1201*Compile_16*Compile_2310)/8.))/4. + \
      (Compile_2315*Compile_2320*Compile_8)/16. - ((Compile_1044 + \
      Compile_2315 + Compile_2323)*Compile_2329*Compile_8)/16. - \
      (Compile_16*Compile_2308*Compile_8*Compile_827)/32.))/2. + \
      (Compile_1133*Compile_348*Compile_809*Compile_817*((Compile_1201*(\
      Compile_2455 - \
      (Compile_1201*Compile_18*Compile_2456)/8.)*Compile_363)/4. + \
      (Compile_2461*Compile_2466*Compile_8)/16. - ((Compile_1046 + \
      Compile_2461 + Compile_2469)*Compile_2475*Compile_8)/16. - \
      (Compile_18*Compile_2453*Compile_8*Compile_827)/32.))/4. + \
      (Compile_1133*Compile_454*Compile_809*Compile_821*((Compile_1201*(\
      Compile_2455 - \
      (Compile_1201*Compile_18*Compile_2524)/8.)*Compile_469)/4. + \
      (Compile_2529*Compile_2534*Compile_8)/16. - ((Compile_1046 + \
      Compile_2529 + Compile_2537)*Compile_2543*Compile_8)/16. - \
      (Compile_18*Compile_2522*Compile_8*Compile_827)/32.))/4. + \
      (Compile_117*Compile_119*Compile_611*Compile_64*Compile_8*Compile_809*\
      Compile_827*(Compile_62*(Compile_788*Compile_833 + \
      2*Compile_577*Compile_838) + 2*Compile_79*(Compile_788*Compile_830 + \
      Compile_577*Compile_841)))/32. + \
      (Compile_117*Compile_119*Compile_611*Compile_64*Compile_8*Compile_809*\
      Compile_827*(Compile_788*(2*Compile_79*Compile_830 + \
      Compile_62*Compile_833) + 2*Compile_577*(Compile_62*Compile_838 + \
      Compile_79*Compile_841)))/32. + \
      (Compile_16*Compile_166*Compile_5*Compile_8*Compile_809*Compile_827*\
      Compile_859)/64. + \
      (Compile_1133*Compile_166*Compile_809*(-0.03125*(Compile_16* Compile_168*Compile_2310*Compile_8) + \
      (Compile_2320*Compile_2338*Compile_8)/16. - \
      (Compile_2329*(Compile_1044 + Compile_2323 + \
      Compile_2338)*Compile_8)/16. + (Compile_1201*(Compile_1376 - \
      (Compile_1201*Compile_16*Compile_2308)/8.)*Compile_827)/4.)*Compile_859)/2. + (Compile_24*Compile_68*(3*lam1H3dUS*pow(Compile_62,4) + \
      Compile_70 + \
      12*Compile_66*Compile_71*Compile_72)*Compile_8*Compile_86)/32. + \
      (Compile_166*Compile_168*Compile_611*Compile_64*Compile_8*Compile_809*\
      Compile_827*(Compile_154*(Compile_788*Compile_864 + \
      2*Compile_577*Compile_868) + 2*Compile_79*(Compile_788*Compile_861 + \
      Compile_577*Compile_871)))/32. + \
      (Compile_166*Compile_168*Compile_611*Compile_64*Compile_8*Compile_809*\
      Compile_827*(Compile_788*(2*Compile_79*Compile_861 + \
      Compile_154*Compile_864) + 2*Compile_577*(Compile_154*Compile_868 + \
      Compile_79*Compile_871)))/32. + \
      (Compile_143*Compile_148*Compile_24*Compile_8*(2*Compile_72*(\
      lam5H3dUS*Compile_127 + Compile_88 + Compile_127*Compile_89) + \
      Compile_154*((Compile_103 + lam1H3dUS*Compile_127)*Compile_154 - \
      2*Compile_126*Compile_72*Compile_89)))/32. + \
      (Compile_117*Compile_119*Compile_166*Compile_168*Compile_64*Compile_8*\
      (-8*m23dUS*lam2H3dUS*lam5H3dUS*Compile_10 + \
      4*m23dUS*lam3H3dUS*lam5H3dUS*Compile_10 - 4*lam5H3dUS*Compile_172 + \
      Compile_173 + Compile_174 + 4*Compile_178*Compile_180 + Compile_183 + \
      Compile_185 + Compile_186 + Compile_187 + \
      4*lam1H3dUS*lam3H3dUS*lam5H3dUS*Compile_191 - \
      4*lam5H3dUS*Compile_190*Compile_191 + Compile_192 + Compile_193 + \
      Compile_195 + Compile_196 + Compile_197 + Compile_198 + Compile_201 + \
      Compile_203 + Compile_204 + Compile_205 + Compile_206 + Compile_209 - \
      8*m12R3dUS*phi1*phi2*Compile_176*Compile_21 + Compile_210 + \
      Compile_211 + Compile_212 + Compile_213 + Compile_214 + Compile_215 + \
      Compile_216 + Compile_226 + Compile_229 + Compile_239 + Compile_240 + \
      Compile_241 + Compile_242 + Compile_243 + Compile_244 + Compile_247 - \
      4*m13dUS*Compile_180*Compile_248 + 8*Compile_176*Compile_34 + \
      4*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_58 - \
      4*lam5H3dUS*Compile_238*Compile_58 + \
      8*m23dUS*lam1H3dUS*lam5H3dUS*Compile_9 - \
      4*m23dUS*lam3H3dUS*lam5H3dUS*Compile_9 + \
      8*lam1H3dUS*lam2H3dUS*lam5H3dUS*Compile_10*Compile_9 - \
      4*lam1H3dUS*lam3H3dUS*lam5H3dUS*Compile_10*Compile_9 - \
      4*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_10*Compile_9 + \
      4*lam1H3dUS*lam4H3dUS*lam5H3dUS*Compile_10*Compile_9 + \
      4*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_10*Compile_9 + \
      2*lam1H3dUS*Compile_10*Compile_217*Compile_9 + \
      2*lam2H3dUS*Compile_10*Compile_217*Compile_9 + \
      2*lam3H3dUS*Compile_10*Compile_217*Compile_9 - \
      2*lam5H3dUS*Compile_10*Compile_217*Compile_9 + \
      2*Compile_10*Compile_221*Compile_9 + \
      2*lam1H3dUS*Compile_10*Compile_231*Compile_9 + \
      2*lam2H3dUS*Compile_10*Compile_231*Compile_9 + \
      2*lam3H3dUS*Compile_10*Compile_231*Compile_9 - \
      10*lam4H3dUS*Compile_10*Compile_231*Compile_9 - \
      6*Compile_10*Compile_236*Compile_9))/2. + \
      (Compile_117*Compile_119*Compile_166*Compile_168*Compile_64*Compile_8*\
      (-4*m23dUS*lam3H3dUS*lam5H3dUS*Compile_10 + Compile_173 + Compile_174 \
      + Compile_183 + Compile_185 + Compile_186 + Compile_187 - \
      4*lam1H3dUS*lam3H3dUS*lam5H3dUS*Compile_191 + \
      lam5H3dUS*Compile_184*Compile_191 + Compile_192 + Compile_193 + \
      Compile_195 + Compile_196 + Compile_197 + Compile_198 + Compile_203 + \
      Compile_204 + Compile_205 + Compile_206 + Compile_209 + Compile_210 + \
      Compile_211 + Compile_212 + Compile_213 + Compile_214 + Compile_215 + \
      Compile_216 + Compile_239 + Compile_240 + Compile_241 + Compile_242 + \
      Compile_243 + Compile_244 - \
      8*m12R3dUS*phi1*phi2*Compile_21*Compile_255 + Compile_256 + \
      Compile_257 + Compile_258 + Compile_259 + Compile_261 + Compile_265 + \
      Compile_270 + Compile_271 + Compile_282 + Compile_284 - \
      4*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_58 + \
      lam5H3dUS*Compile_184*Compile_58 - 4*m13dUS*Compile_248*Compile_71 + \
      4*m23dUS*lam3H3dUS*lam5H3dUS*Compile_9 + \
      4*lam1H3dUS*lam3H3dUS*lam5H3dUS*Compile_10*Compile_9 + \
      4*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_10*Compile_9 + \
      12*lam1H3dUS*lam4H3dUS*lam5H3dUS*Compile_10*Compile_9 + \
      12*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_10*Compile_9 - \
      4*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_10*Compile_9 - \
      2*lam5H3dUS*Compile_10*Compile_184*Compile_9 + \
      6*lam1H3dUS*Compile_10*Compile_217*Compile_9 + \
      6*lam2H3dUS*Compile_10*Compile_217*Compile_9 - \
      2*lam3H3dUS*Compile_10*Compile_217*Compile_9 - \
      6*lam5H3dUS*Compile_10*Compile_217*Compile_9 + \
      6*lam1H3dUS*Compile_10*Compile_231*Compile_9 + \
      6*lam2H3dUS*Compile_10*Compile_231*Compile_9 - \
      2*lam3H3dUS*Compile_10*Compile_231*Compile_9 - \
      2*Compile_10*Compile_236*Compile_9))/2. + \
      (Compile_325*Compile_348*Compile_363*Compile_611*Compile_8*Compile_809*Compile_827*(Compile_788*(8*Compile_291*Compile_889 + \
      Compile_369*Compile_892) + 8*Compile_577*(Compile_369*Compile_897 + \
      Compile_291*Compile_900)))/64. + \
      (Compile_325*Compile_348*Compile_363*Compile_611*Compile_8*Compile_809*Compile_827*(Compile_369*(Compile_788*Compile_892 + \
      8*Compile_577*Compile_897) + 8*Compile_291*(Compile_788*Compile_889 + \
      Compile_577*Compile_900)))/64. + \
      (Compile_325*Compile_454*Compile_469*Compile_611*Compile_8*Compile_809*Compile_827*(Compile_788*(8*Compile_291*Compile_917 + \
      Compile_431*Compile_920) + 8*Compile_577*(Compile_431*Compile_924 + \
      Compile_291*Compile_927)))/64. + \
      (Compile_325*Compile_454*Compile_469*Compile_611*Compile_8*Compile_809*Compile_827*(Compile_431*(Compile_788*Compile_920 + \
      8*Compile_577*Compile_924) + 8*Compile_291*(Compile_788*Compile_917 + \
      Compile_577*Compile_927)))/64. + \
      Compile_325*Compile_348*Compile_363*Compile_454*Compile_469*Compile_8*\
      (-8*m23dUS*lam3H3dUS*lam4H3dUS*Compile_10 + Compile_173 + Compile_174 \
      + Compile_183 + Compile_185 + Compile_186 - \
      8*lam1H3dUS*lam3H3dUS*lam4H3dUS*Compile_191 - \
      2*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_191 + \
      3*lam4H3dUS*Compile_184*Compile_191 + Compile_192 + Compile_193 + \
      Compile_195 + Compile_196 + Compile_201 + Compile_203 + Compile_204 + \
      Compile_205 + Compile_209 + Compile_210 + Compile_211 + Compile_212 + \
      Compile_213 - 4*m23dUS*Compile_10*Compile_217 - \
      4*lam1H3dUS*Compile_191*Compile_217 + \
      3*lam3H3dUS*Compile_191*Compile_217 - \
      lam5H3dUS*Compile_191*Compile_217 + Compile_191*Compile_221 + \
      Compile_226 + Compile_229 + 4*m23dUS*Compile_10*Compile_231 + \
      4*lam1H3dUS*Compile_191*Compile_231 - \
      lam3H3dUS*Compile_191*Compile_231 - lam4H3dUS*Compile_191*Compile_231 \
      + Compile_191*Compile_236 + Compile_239 + Compile_240 + Compile_241 + \
      Compile_242 + Compile_247 - \
      16*m12R3dUS*lam5H3dUS*phi1*phi2*Compile_255 + Compile_256 + \
      Compile_257 + Compile_258 + Compile_259 + Compile_261 + Compile_265 + \
      Compile_270 + Compile_271 + Compile_282 + Compile_284 - \
      8*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_58 - \
      2*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_58 + \
      3*lam4H3dUS*Compile_184*Compile_58 - \
      4*lam2H3dUS*Compile_217*Compile_58 + \
      3*lam3H3dUS*Compile_217*Compile_58 - lam5H3dUS*Compile_217*Compile_58 \
      + Compile_221*Compile_58 + 4*lam2H3dUS*Compile_231*Compile_58 - \
      lam3H3dUS*Compile_231*Compile_58 - lam4H3dUS*Compile_231*Compile_58 + \
      Compile_236*Compile_58 + 8*m23dUS*lam3H3dUS*lam4H3dUS*Compile_9 + \
      8*lam1H3dUS*lam3H3dUS*lam4H3dUS*Compile_10*Compile_9 + \
      8*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_10*Compile_9 - \
      6*lam4H3dUS*Compile_10*Compile_184*Compile_9 + \
      4*m23dUS*Compile_217*Compile_9 + \
      4*lam1H3dUS*Compile_10*Compile_217*Compile_9 + \
      4*lam2H3dUS*Compile_10*Compile_217*Compile_9 - \
      6*lam3H3dUS*Compile_10*Compile_217*Compile_9 + \
      2*lam5H3dUS*Compile_10*Compile_217*Compile_9 - \
      4*m23dUS*Compile_231*Compile_9 + \
      20*lam1H3dUS*Compile_10*Compile_231*Compile_9 + \
      20*lam2H3dUS*Compile_10*Compile_231*Compile_9 - \
      6*lam3H3dUS*Compile_10*Compile_231*Compile_9 - \
      10*Compile_10*Compile_236*Compile_9 - \
      4*m13dUS*Compile_71*(Compile_568 + Compile_569 + Compile_570 + \
      Compile_571 + Compile_95 + Compile_96 + Compile_97 + Compile_98 + \
      Compile_99)) + \
      Compile_611*Compile_637*Compile_655*Compile_8*Compile_809*Compile_827*\
      (Compile_256 - 16*m12R3dUS*phi1*phi2*Compile_255*Compile_71 + \
      Compile_71*(4*Compile_172 + 4*Compile_178 - \
      12*lam1H3dUS*lam3H3dUS*Compile_191 - \
      12*lam1H3dUS*lam4H3dUS*Compile_191 + \
      2*lam3H3dUS*lam4H3dUS*Compile_191 - \
      12*lam1H3dUS*lam5H3dUS*Compile_191 + \
      2*lam3H3dUS*lam5H3dUS*Compile_191 + 2*lam4H3dUS*lam5H3dUS*Compile_191 \
      + Compile_184*Compile_191 + 36*Compile_190*Compile_191 + \
      Compile_191*Compile_217 + Compile_191*Compile_231 - \
      12*lam2H3dUS*lam3H3dUS*Compile_58 - 12*lam2H3dUS*lam4H3dUS*Compile_58 \
      + 2*lam3H3dUS*lam4H3dUS*Compile_58 - \
      12*lam2H3dUS*lam5H3dUS*Compile_58 + 2*lam3H3dUS*lam5H3dUS*Compile_58 \
      + 2*lam4H3dUS*lam5H3dUS*Compile_58 + Compile_184*Compile_58 + \
      Compile_217*Compile_58 + Compile_231*Compile_58 + \
      36*Compile_238*Compile_58 - \
      72*lam1H3dUS*lam2H3dUS*Compile_10*Compile_9 + \
      36*lam1H3dUS*lam3H3dUS*Compile_10*Compile_9 + \
      36*lam2H3dUS*lam3H3dUS*Compile_10*Compile_9 + \
      36*lam1H3dUS*lam4H3dUS*Compile_10*Compile_9 + \
      36*lam2H3dUS*lam4H3dUS*Compile_10*Compile_9 - \
      20*lam3H3dUS*lam4H3dUS*Compile_10*Compile_9 + \
      36*lam1H3dUS*lam5H3dUS*Compile_10*Compile_9 + \
      36*lam2H3dUS*lam5H3dUS*Compile_10*Compile_9 - \
      20*lam3H3dUS*lam5H3dUS*Compile_10*Compile_9 - \
      20*lam4H3dUS*lam5H3dUS*Compile_10*Compile_9 - \
      10*Compile_10*Compile_184*Compile_9 - \
      10*Compile_10*Compile_217*Compile_9 - \
      10*Compile_10*Compile_231*Compile_9 - 4*m23dUS*(Compile_31 + \
      Compile_33 + Compile_569 + Compile_571 + Compile_657 + Compile_658 - \
      lam4H3dUS*Compile_9 + Compile_982) - 4*m13dUS*(Compile_568 + \
      Compile_570 + Compile_580 + Compile_583 + Compile_584 + Compile_585 + \
      Compile_95 + Compile_97 + Compile_99))) + \
      (Compile_24*Compile_68*Compile_8*Compile_86*(2*Compile_72*(lam5H3dUS*\
      Compile_66 + Compile_88 + Compile_66*Compile_89) + \
      Compile_62*(Compile_62*(Compile_103 + lam1H3dUS*Compile_66) - \
      2*Compile_72*Compile_89*(Compile_100 + Compile_94 + Compile_95 + \
      Compile_96 + Compile_97 + Compile_98 + Compile_99))))/32. + \
      (Compile_11*Compile_13*Compile_7*Compile_8*Compile_992)/128. + \
      (Compile_11*Compile_2*Compile_3*Compile_7*Compile_8*Compile_992)/128. \
      + (Compile_1011*Compile_1132*Compile_1165*Compile_1242*Compile_166*(\
      Compile_1166 + Compile_1259 + Compile_1260 + Compile_1263 + \
      Compile_1265 + Compile_1268 + Compile_1276 - \
      (Compile_1283*(Compile_1044 + Compile_1045 + Compile_1270 - \
      (Compile_11*Compile_143*Compile_166*Compile_5)/2.)*Compile_8)/16.)*\
      Compile_996)/2. + \
      (Compile_1011*Compile_1132*Compile_1165*Compile_1242*Compile_166*(\
      Compile_1166 + Compile_1259 + Compile_1260 + Compile_1263 + \
      Compile_1265 + Compile_1268 + Compile_1276 - \
      (Compile_1283*(Compile_1045 + Compile_1046 + Compile_1267 - \
      (Compile_11*Compile_143*Compile_166*Compile_6)/2.)*Compile_8)/16.)*\
      Compile_996)/2. + \
      (Compile_1011*Compile_1132*Compile_1146*Compile_1165*Compile_117*(\
      Compile_1166 + Compile_1167 + Compile_1168 + Compile_1171 + \
      Compile_1173 + Compile_1176 + Compile_1184 - \
      (Compile_1191*(Compile_1044 + Compile_1045 + Compile_1178 - \
      (Compile_11*Compile_117*Compile_5*Compile_68)/2.)*Compile_8)/16.)*\
      Compile_996)/2. + \
      (Compile_1011*Compile_1132*Compile_1146*Compile_1165*Compile_117*(\
      Compile_1166 + Compile_1167 + Compile_1168 + Compile_1171 + \
      Compile_1173 + Compile_1176 + Compile_1184 - \
      (Compile_1191*(Compile_1045 + Compile_1046 + Compile_1175 - \
      (Compile_11*Compile_117*Compile_6*Compile_68)/2.)*Compile_8)/16.)*\
      Compile_996)/2. - (Compile_1498*Compile_632*pow(Compile_1520 + \
      Compile_1512*Compile_661,2)*Compile_8*(0.5 + \
      log((mu3US*Compile_1523)/3.)))/48. - (Compile_1498*pow(Compile_2178 + \
      Compile_2171*Compile_788,2)*Compile_8*Compile_803*(0.5 + \
      log((mu3US*Compile_2181)/3.)))/48. - \
      (Compile_1132*Compile_117*Compile_1389*Compile_166*pow(Compile_1390 + \
      Compile_1391 + Compile_1392 + Compile_1393 + Compile_1394 + \
      Compile_1395 + Compile_1396 + Compile_1398 + Compile_1399 + \
      Compile_1400 + Compile_1401 - \
      phi2*Compile_322,2)*Compile_325*Compile_348*Compile_60*Compile_8*(0.5 \
      + log(mu3US/(Compile_119 + Compile_168 + Compile_363))))/4. - \
      (Compile_1132*Compile_117*Compile_1389*Compile_166*pow(Compile_1390 + \
      Compile_1391 + Compile_1392 + Compile_1393 + Compile_1394 + \
      Compile_1395 + Compile_1396 + Compile_1398 + Compile_1399 + \
      Compile_1400 + Compile_1401 + \
      phi2*Compile_322,2)*Compile_325*Compile_454*Compile_60*Compile_8*(0.5 \
      + log(mu3US/(Compile_119 + Compile_168 + Compile_469))))/4. + \
      (Compile_1537*Compile_611*Compile_637*Compile_996*(Compile_1538 + \
      Compile_1549 - (Compile_11*(Compile_1540 + \
      Compile_1541)*Compile_5*Compile_8)/64. + \
      (Compile_11*Compile_16*Compile_5*Compile_655*Compile_8)/64. + \
      (Compile_1551*Compile_8*(0.5 + log(mu3US/(Compile_1051 + \
      Compile_655))))/8. - ((Compile_1551 + Compile_1558 - \
      (Compile_11*Compile_5*Compile_619*Compile_637)/2.)*Compile_8*(0.5 + \
      log(mu3US/(Compile_16 + Compile_655))))/16.))/2. - \
      (Compile_1132*Compile_117*Compile_166*Compile_611*Compile_637*Compile_8*pow(-16*m12R3dUS*m13dUS*lam1H3dUS*phi1 + \
      16*m12R3dUS*m23dUS*lam1H3dUS*phi1 + 8*m12R3dUS*m13dUS*lam3H3dUS*phi1 \
      - 8*m12R3dUS*m23dUS*lam3H3dUS*phi1 + 8*m12R3dUS*m13dUS*lam4H3dUS*phi1 \
      - 8*m12R3dUS*m23dUS*lam4H3dUS*phi1 + 8*m12R3dUS*m13dUS*lam5H3dUS*phi1 \
      - 8*m12R3dUS*m23dUS*lam5H3dUS*phi1 + 8*m13dUS*m23dUS*lam4H3dUS*phi2 + \
      8*m13dUS*m23dUS*lam5H3dUS*phi2 + \
      48*m12R3dUS*lam1H3dUS*lam2H3dUS*phi1*Compile_10 - \
      8*m12R3dUS*lam1H3dUS*lam3H3dUS*phi1*Compile_10 + \
      8*m12R3dUS*lam2H3dUS*lam3H3dUS*phi1*Compile_10 - \
      8*m12R3dUS*lam1H3dUS*lam4H3dUS*phi1*Compile_10 + \
      40*m12R3dUS*lam2H3dUS*lam4H3dUS*phi1*Compile_10 - \
      16*m12R3dUS*lam3H3dUS*lam4H3dUS*phi1*Compile_10 - \
      8*m12R3dUS*lam1H3dUS*lam5H3dUS*phi1*Compile_10 + \
      40*m12R3dUS*lam2H3dUS*lam5H3dUS*phi1*Compile_10 - \
      16*m12R3dUS*lam3H3dUS*lam5H3dUS*phi1*Compile_10 + \
      32*m12R3dUS*lam1H3dUS*lam3H3dUS*Compile_1136 + \
      16*m12R3dUS*lam1H3dUS*lam4H3dUS*Compile_1136 - \
      8*m12R3dUS*lam3H3dUS*lam4H3dUS*Compile_1136 + \
      16*m12R3dUS*lam1H3dUS*lam5H3dUS*Compile_1136 - \
      8*m12R3dUS*lam3H3dUS*lam5H3dUS*Compile_1136 + \
      16*m13dUS*lam2H3dUS*lam4H3dUS*Compile_1397 - \
      16*m23dUS*lam2H3dUS*lam4H3dUS*Compile_1397 - \
      4*m13dUS*lam3H3dUS*lam4H3dUS*Compile_1397 + \
      4*m23dUS*lam3H3dUS*lam4H3dUS*Compile_1397 + \
      16*m13dUS*lam2H3dUS*lam5H3dUS*Compile_1397 - \
      16*m23dUS*lam2H3dUS*lam5H3dUS*Compile_1397 - \
      4*m13dUS*lam3H3dUS*lam5H3dUS*Compile_1397 + \
      4*m23dUS*lam3H3dUS*lam5H3dUS*Compile_1397 - \
      4*m13dUS*lam4H3dUS*lam5H3dUS*Compile_1397 + \
      4*m23dUS*lam4H3dUS*lam5H3dUS*Compile_1397 - \
      4*lam4H3dUS*phi2*Compile_172 - 4*lam5H3dUS*phi2*Compile_172 - \
      4*lam4H3dUS*phi2*Compile_178 - 4*lam5H3dUS*phi2*Compile_178 + \
      Compile_1804 + 8*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_1810 + \
      8*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_1810 + \
      4*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_1810 - \
      2*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_1810 + Compile_1823 + \
      Compile_1824 + Compile_1825 + Compile_1826 + Compile_1827 + \
      Compile_1828 + Compile_1829 + Compile_1830 + Compile_1831 + \
      Compile_1832 + Compile_1833 + Compile_1834 + Compile_1835 + \
      Compile_1836 - 12*m12R3dUS*phi1*Compile_10*Compile_184 - \
      4*m12R3dUS*Compile_1136*Compile_184 - \
      lam4H3dUS*Compile_1810*Compile_184 - \
      lam5H3dUS*Compile_1810*Compile_184 - \
      48*m12R3dUS*Compile_1136*Compile_190 - \
      16*lam1H3dUS*lam3H3dUS*lam4H3dUS*phi2*Compile_191 - \
      16*lam1H3dUS*lam3H3dUS*lam5H3dUS*phi2*Compile_191 - \
      20*lam1H3dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_191 + \
      10*lam3H3dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_191 + \
      5*lam4H3dUS*phi2*Compile_184*Compile_191 + \
      5*lam5H3dUS*phi2*Compile_184*Compile_191 + \
      12*lam4H3dUS*phi2*Compile_190*Compile_191 + \
      12*lam5H3dUS*phi2*Compile_190*Compile_191 - \
      2*m13dUS*Compile_1397*Compile_217 + 2*m23dUS*Compile_1397*Compile_217 \
      + 2*lam2H3dUS*Compile_1810*Compile_217 - \
      lam3H3dUS*Compile_1810*Compile_217 - \
      10*lam1H3dUS*phi2*Compile_191*Compile_217 + \
      5*lam3H3dUS*phi2*Compile_191*Compile_217 - \
      2*m13dUS*Compile_1397*Compile_231 + 2*m23dUS*Compile_1397*Compile_231 \
      + 2*lam2H3dUS*Compile_1810*Compile_231 - \
      lam3H3dUS*Compile_1810*Compile_231 - \
      10*lam1H3dUS*phi2*Compile_191*Compile_231 + \
      5*lam3H3dUS*phi2*Compile_191*Compile_231 - \
      12*lam4H3dUS*Compile_1810*Compile_238 - \
      12*lam5H3dUS*Compile_1810*Compile_238 - 32*lam2H3dUS*phi2*Compile_34 \
      + 16*lam3H3dUS*phi2*Compile_34 - \
      8*m13dUS*lam1H3dUS*lam4H3dUS*phi2*Compile_9 + \
      8*m23dUS*lam1H3dUS*lam4H3dUS*phi2*Compile_9 - \
      8*m13dUS*lam3H3dUS*lam4H3dUS*phi2*Compile_9 + \
      8*m23dUS*lam3H3dUS*lam4H3dUS*phi2*Compile_9 - \
      8*m13dUS*lam1H3dUS*lam5H3dUS*phi2*Compile_9 + \
      8*m23dUS*lam1H3dUS*lam5H3dUS*phi2*Compile_9 - \
      8*m13dUS*lam3H3dUS*lam5H3dUS*phi2*Compile_9 + \
      8*m23dUS*lam3H3dUS*lam5H3dUS*phi2*Compile_9 - \
      12*m13dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_9 + \
      12*m23dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_9 - \
      4*lam1H3dUS*lam3H3dUS*lam4H3dUS*Compile_1397*Compile_9 - \
      4*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_1397*Compile_9 - \
      4*lam1H3dUS*lam3H3dUS*lam5H3dUS*Compile_1397*Compile_9 - \
      4*lam2H3dUS*lam3H3dUS*lam5H3dUS*Compile_1397*Compile_9 - \
      20*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_1397*Compile_9 + \
      8*lam3H3dUS*lam4H3dUS*lam5H3dUS*Compile_1397*Compile_9 + \
      4*lam4H3dUS*Compile_1397*Compile_184*Compile_9 + \
      4*lam5H3dUS*Compile_1397*Compile_184*Compile_9 - \
      6*m13dUS*phi2*Compile_217*Compile_9 + \
      6*m23dUS*phi2*Compile_217*Compile_9 + \
      2*lam1H3dUS*Compile_1397*Compile_217*Compile_9 - \
      10*lam2H3dUS*Compile_1397*Compile_217*Compile_9 + \
      4*lam3H3dUS*Compile_1397*Compile_217*Compile_9 - \
      6*m13dUS*phi2*Compile_231*Compile_9 + \
      6*m23dUS*phi2*Compile_231*Compile_9 + \
      2*lam1H3dUS*Compile_1397*Compile_231*Compile_9 - \
      10*lam2H3dUS*Compile_1397*Compile_231*Compile_9 + \
      4*lam3H3dUS*Compile_1397*Compile_231*Compile_9,2)*(0.5 + \
      log(mu3US/(Compile_119 + Compile_168 + Compile_655))))/4. + \
      (Compile_1537*Compile_611*Compile_637*Compile_996*(Compile_1549 + \
      Compile_1569 - (Compile_11*(Compile_1541 + \
      Compile_1571)*Compile_6*Compile_8)/64. + \
      (Compile_11*Compile_18*Compile_6*Compile_655*Compile_8)/64. + \
      (Compile_1575*Compile_8*(0.5 + log(mu3US/(Compile_1052 + \
      Compile_655))))/8. - ((Compile_1575 + Compile_1582 - \
      (Compile_11*Compile_6*Compile_619*Compile_637)/2.)*Compile_8*(0.5 + \
      log(mu3US/(Compile_18 + Compile_655))))/16.))/4. - \
      2*Compile_325*Compile_348*Compile_454*Compile_611*Compile_637*Compile_8*pow(Compile_1804 + Compile_2020 + Compile_2021 + Compile_2022 + \
      Compile_2023 + Compile_2024 + Compile_2025 + Compile_2026 + \
      Compile_2027 + Compile_2028 + Compile_2029 + Compile_2030 + \
      Compile_2031 + Compile_2032 + Compile_2033 + Compile_2034 + \
      Compile_2035 + Compile_2036 + Compile_2037 + Compile_2038 + \
      Compile_2039 + Compile_2040 + Compile_2041 + Compile_2042 + \
      Compile_2043 + Compile_2044 + Compile_2045 + Compile_2046 + \
      Compile_2047 + Compile_2048 + Compile_2049 + Compile_2050 + \
      Compile_2051 + Compile_2052 + Compile_2053 + Compile_2054 + \
      Compile_2055 + Compile_2056 + Compile_2057 + Compile_2058 + \
      Compile_2059 + Compile_2060 + Compile_2061 + Compile_2062 + \
      Compile_2063 + Compile_2064 + Compile_2065 + Compile_2066 + \
      Compile_2067 + Compile_2068 + Compile_2069 + Compile_2070 + \
      Compile_2071 + Compile_2072 + Compile_2073 + Compile_2074 + \
      Compile_2075 + Compile_2076 + Compile_2077 + Compile_2078 + \
      Compile_2079 + Compile_2080 + Compile_2081 + Compile_2082 + \
      Compile_2083 + Compile_2084 + Compile_2085 + Compile_2086 + \
      Compile_2087 + Compile_2088 + Compile_2089 + Compile_2090 + \
      Compile_2091 + Compile_2092 + Compile_2093 + Compile_2094 + \
      Compile_2095 - 4*m12R3dUS*lam1H3dUS*phi1*Compile_608 + \
      2*m12R3dUS*lam3H3dUS*phi1*Compile_608 + \
      2*m12R3dUS*lam4H3dUS*phi1*Compile_608 - \
      2*m12R3dUS*lam5H3dUS*phi1*Compile_608 - \
      2*m13dUS*lam5H3dUS*phi2*Compile_608 + \
      2*m23dUS*lam5H3dUS*phi2*Compile_608 + \
      2*lam2H3dUS*lam5H3dUS*Compile_1397*Compile_608 - \
      lam3H3dUS*lam5H3dUS*Compile_1397*Compile_608 - \
      lam4H3dUS*lam5H3dUS*Compile_1397*Compile_608 + \
      Compile_1397*Compile_231*Compile_608 + \
      2*lam1H3dUS*lam5H3dUS*phi2*Compile_608*Compile_9 - \
      lam3H3dUS*lam5H3dUS*phi2*Compile_608*Compile_9 - \
      lam4H3dUS*lam5H3dUS*phi2*Compile_608*Compile_9 + \
      phi2*Compile_231*Compile_608*Compile_9,2)*(0.5 + \
      log(mu3US/(Compile_363 + Compile_469 + Compile_655))) + \
      (Compile_2188*Compile_611*Compile_809*Compile_996*(Compile_1538 + \
      Compile_2198 - (Compile_11*(Compile_1540 + \
      Compile_2190)*Compile_5*Compile_8)/64. + \
      (Compile_11*Compile_16*Compile_5*Compile_8*Compile_827)/64. + \
      (Compile_2200*Compile_8*(0.5 + log(mu3US/(Compile_1051 + \
      Compile_827))))/8. - (Compile_8*(Compile_1558 + Compile_2200 - \
      (Compile_11*Compile_5*Compile_792*Compile_809)/2.)*(0.5 + \
      log(mu3US/(Compile_16 + Compile_827))))/16.))/2. - \
      (Compile_1132*Compile_117*Compile_166*Compile_611*Compile_8*Compile_809*pow(16*m12R3dUS*m13dUS*lam1H3dUS*phi1 - \
      16*m12R3dUS*m23dUS*lam1H3dUS*phi1 - 8*m12R3dUS*m13dUS*lam3H3dUS*phi1 \
      + 8*m12R3dUS*m23dUS*lam3H3dUS*phi1 - 8*m12R3dUS*m13dUS*lam4H3dUS*phi1 \
      + 8*m12R3dUS*m23dUS*lam4H3dUS*phi1 - 8*m12R3dUS*m13dUS*lam5H3dUS*phi1 \
      + 8*m12R3dUS*m23dUS*lam5H3dUS*phi1 - 8*m13dUS*m23dUS*lam4H3dUS*phi2 - \
      48*m12R3dUS*lam1H3dUS*lam2H3dUS*phi1*Compile_10 + \
      8*m12R3dUS*lam1H3dUS*lam3H3dUS*phi1*Compile_10 - \
      8*m12R3dUS*lam2H3dUS*lam3H3dUS*phi1*Compile_10 + \
      8*m12R3dUS*lam1H3dUS*lam4H3dUS*phi1*Compile_10 - \
      40*m12R3dUS*lam2H3dUS*lam4H3dUS*phi1*Compile_10 + \
      16*m12R3dUS*lam3H3dUS*lam4H3dUS*phi1*Compile_10 + \
      8*m12R3dUS*lam1H3dUS*lam5H3dUS*phi1*Compile_10 - \
      40*m12R3dUS*lam2H3dUS*lam5H3dUS*phi1*Compile_10 + \
      16*m12R3dUS*lam3H3dUS*lam5H3dUS*phi1*Compile_10 - \
      32*m12R3dUS*lam1H3dUS*lam3H3dUS*Compile_1136 + \
      8*m12R3dUS*lam3H3dUS*lam4H3dUS*Compile_1136 - \
      16*m12R3dUS*lam1H3dUS*lam5H3dUS*Compile_1136 + \
      8*m12R3dUS*lam3H3dUS*lam5H3dUS*Compile_1136 - \
      16*m13dUS*lam2H3dUS*lam4H3dUS*Compile_1397 + \
      16*m23dUS*lam2H3dUS*lam4H3dUS*Compile_1397 + \
      4*m13dUS*lam3H3dUS*lam4H3dUS*Compile_1397 - \
      4*m23dUS*lam3H3dUS*lam4H3dUS*Compile_1397 + \
      4*lam4H3dUS*phi2*Compile_172 + 4*lam4H3dUS*phi2*Compile_178 - \
      8*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_1810 - \
      4*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_1810 + Compile_1823 + \
      Compile_1824 + Compile_1825 + Compile_1826 + Compile_1827 + \
      Compile_1828 + Compile_1829 + Compile_1830 + Compile_1831 + \
      Compile_1832 + Compile_1833 + Compile_1834 + Compile_1835 + \
      Compile_1836 + 12*m12R3dUS*phi1*Compile_10*Compile_184 + \
      4*m12R3dUS*Compile_1136*Compile_184 + \
      lam4H3dUS*Compile_1810*Compile_184 + \
      48*m12R3dUS*Compile_1136*Compile_190 + \
      16*lam1H3dUS*lam3H3dUS*lam4H3dUS*phi2*Compile_191 + \
      20*lam1H3dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_191 - \
      5*lam4H3dUS*phi2*Compile_184*Compile_191 - \
      12*lam4H3dUS*phi2*Compile_190*Compile_191 + Compile_2031 + \
      Compile_2041 + Compile_2042 + Compile_2043 + Compile_2044 + \
      Compile_2045 + Compile_2046 + Compile_2047 + Compile_2052 + \
      Compile_2053 + Compile_2054 + Compile_2056 + Compile_2073 + \
      Compile_2074 + Compile_2075 + Compile_2076 + Compile_2077 + \
      Compile_2078 + Compile_2079 + Compile_2080 + Compile_2081 + \
      Compile_2083 + Compile_2088 + Compile_2089 + Compile_2090 + \
      Compile_2092 + 2*m13dUS*Compile_1397*Compile_217 - \
      2*m23dUS*Compile_1397*Compile_217 - \
      2*lam2H3dUS*Compile_1810*Compile_217 + \
      lam3H3dUS*Compile_1810*Compile_217 + \
      10*lam1H3dUS*phi2*Compile_191*Compile_217 - \
      5*lam3H3dUS*phi2*Compile_191*Compile_217 + \
      2*m13dUS*Compile_1397*Compile_231 - 2*m23dUS*Compile_1397*Compile_231 \
      - 2*lam2H3dUS*Compile_1810*Compile_231 + \
      lam3H3dUS*Compile_1810*Compile_231 + \
      10*lam1H3dUS*phi2*Compile_191*Compile_231 - \
      5*lam3H3dUS*phi2*Compile_191*Compile_231 + \
      12*lam4H3dUS*Compile_1810*Compile_238 + 32*lam2H3dUS*phi2*Compile_34 \
      - 16*lam3H3dUS*phi2*Compile_34 + \
      8*m13dUS*lam1H3dUS*lam4H3dUS*phi2*Compile_9 - \
      8*m23dUS*lam1H3dUS*lam4H3dUS*phi2*Compile_9 + \
      8*m13dUS*lam3H3dUS*lam4H3dUS*phi2*Compile_9 - \
      8*m23dUS*lam3H3dUS*lam4H3dUS*phi2*Compile_9 + \
      12*m13dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_9 - \
      12*m23dUS*lam4H3dUS*lam5H3dUS*phi2*Compile_9 + \
      4*lam1H3dUS*lam3H3dUS*lam4H3dUS*Compile_1397*Compile_9 + \
      4*lam2H3dUS*lam3H3dUS*lam4H3dUS*Compile_1397*Compile_9 - \
      4*lam1H3dUS*lam4H3dUS*lam5H3dUS*Compile_1397*Compile_9 + \
      20*lam2H3dUS*lam4H3dUS*lam5H3dUS*Compile_1397*Compile_9 - \
      4*lam4H3dUS*Compile_1397*Compile_184*Compile_9 + \
      6*m13dUS*phi2*Compile_217*Compile_9 - \
      6*m23dUS*phi2*Compile_217*Compile_9 - \
      2*lam1H3dUS*Compile_1397*Compile_217*Compile_9 + \
      10*lam2H3dUS*Compile_1397*Compile_217*Compile_9 - \
      4*lam3H3dUS*Compile_1397*Compile_217*Compile_9 + \
      6*m13dUS*phi2*Compile_231*Compile_9 - \
      6*m23dUS*phi2*Compile_231*Compile_9 - \
      2*lam1H3dUS*Compile_1397*Compile_231*Compile_9 + \
      10*lam2H3dUS*Compile_1397*Compile_231*Compile_9 - \
      4*lam3H3dUS*Compile_1397*Compile_231*Compile_9,2)*(0.5 + \
      log(mu3US/(Compile_119 + Compile_168 + Compile_827))))/4. + \
      (Compile_2188*Compile_611*Compile_809*Compile_996*(Compile_1569 + \
      Compile_2198 - (Compile_11*(Compile_1571 + \
      Compile_2190)*Compile_6*Compile_8)/64. + \
      (Compile_11*Compile_18*Compile_6*Compile_8*Compile_827)/64. + \
      (Compile_2221*Compile_8*(0.5 + log(mu3US/(Compile_1052 + \
      Compile_827))))/8. - (Compile_8*(Compile_1582 + Compile_2221 - \
      (Compile_11*Compile_6*Compile_792*Compile_809)/2.)*(0.5 + \
      log(mu3US/(Compile_18 + Compile_827))))/16.))/4. - \
      2*Compile_325*Compile_348*Compile_454*Compile_611*Compile_8*Compile_809*pow(Compile_1804 + Compile_1827 + Compile_1828 + Compile_1831 + \
      Compile_1832 + Compile_1835 + Compile_1836 + Compile_2020 + \
      Compile_2021 + Compile_2022 + Compile_2023 + Compile_2024 + \
      Compile_2025 + Compile_2026 + Compile_2027 + Compile_2028 + \
      Compile_2029 + Compile_2030 + Compile_2031 + Compile_2032 + \
      Compile_2033 + Compile_2034 + Compile_2035 + Compile_2036 + \
      Compile_2037 + Compile_2038 + Compile_2039 + Compile_2040 + \
      Compile_2041 + Compile_2042 + Compile_2043 + Compile_2044 + \
      Compile_2045 + Compile_2046 + Compile_2047 + Compile_2048 + \
      Compile_2049 + Compile_2050 + Compile_2051 + Compile_2052 + \
      Compile_2053 + Compile_2054 + Compile_2055 + Compile_2056 + \
      Compile_2057 + Compile_2058 + Compile_2059 + Compile_2060 + \
      Compile_2061 + Compile_2062 + Compile_2063 + Compile_2064 + \
      Compile_2065 + Compile_2066 + Compile_2067 + Compile_2068 + \
      Compile_2069 + Compile_2070 + Compile_2071 + Compile_2072 + \
      Compile_2073 + Compile_2074 + Compile_2075 + Compile_2076 + \
      Compile_2077 + Compile_2078 + Compile_2079 + Compile_2080 + \
      Compile_2081 + Compile_2082 + Compile_2083 + Compile_2084 + \
      Compile_2085 + Compile_2086 + Compile_2087 + Compile_2088 + \
      Compile_2089 + Compile_2090 + Compile_2091 + Compile_2092 + \
      Compile_2093 + Compile_2094 + Compile_2095 + \
      4*m12R3dUS*lam1H3dUS*phi1*Compile_608 - \
      2*m12R3dUS*lam3H3dUS*phi1*Compile_608 - \
      2*m12R3dUS*lam4H3dUS*phi1*Compile_608 + \
      2*m12R3dUS*lam5H3dUS*phi1*Compile_608 + \
      lam4H3dUS*lam5H3dUS*Compile_1397*Compile_608 - \
      Compile_1397*Compile_231*Compile_608 + \
      lam4H3dUS*lam5H3dUS*phi2*Compile_608*Compile_9 - \
      phi2*Compile_231*Compile_608*Compile_9,2)*(0.5 + \
      log(mu3US/(Compile_363 + Compile_469 + Compile_827))) + \
      (4*Compile_2*Compile_7*Compile_996*((-89*Compile_13*Compile_8*Compile_997)/2048. + (3*Compile_13*Compile_8*Compile_992*Compile_997)/256. + \
      (3*Compile_13*Compile_8*Compile_997*(0.5 + \
      log(mu3US*Compile_989)))/64.))/(3.*pow(g23dUS,2));
    }


  };