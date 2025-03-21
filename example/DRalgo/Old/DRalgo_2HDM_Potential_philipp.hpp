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

      return 0;
    }

    std::complex<double> V1( Eigen::VectorXd phi, double T, std::vector<std::complex<double>> par ) {

      return 0;
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

      std::complex<double> Compile_12 = pow(m13dUS,2);
      std::complex<double> Compile_13 = pow(m23dUS,2);
      std::complex<double> Compile_15 = pow(m12R3dUS,4);
      std::complex<double> Compile_16 = 4.*Compile_15;
      std::complex<double> Compile_17 = -1.*Compile_13;
      std::complex<double> Compile_18 = Compile_12 + Compile_17;
      std::complex<double> Compile_20 = pow(Compile_18,2);
      std::complex<double> Compile_22 = Compile_16 + Compile_20;
      std::complex<double> Compile_24 = sqrt(Compile_22);
      std::complex<double> Compile_44 = pow(phi1,2);
      std::complex<double> Compile_53 = pow(phi2,2);
      std::complex<double> Compile_61 = Compile_44 + Compile_53;
      std::complex<double> Compile_60 = pow(g23dUS,4);
      std::complex<double> Compile_64 = pow(g23dUS,2);
      std::complex<double> Compile_73 = 1./sqrt(Compile_61);
      std::complex<double> Compile_72 = 1./g23dUS;
      std::complex<double> Compile_74 = 2.*mu3US*Compile_72*Compile_73;
      std::complex<double> Compile_75 = log(Compile_74);
      std::complex<double> Compile_79 = 1./sqrt(Compile_64);
      std::complex<double> Compile_80 = 2.*mu3US*Compile_73*Compile_79;
      std::complex<double> Compile_129 = log(Compile_80);
      std::complex<double> Compile_260 = pow(g23dUS,6);
      std::complex<double> Compile_269 = sqrt(Compile_61);
      std::complex<double> Compile_51 = pow(m12R3dUS,2);
      std::complex<double> Compile_337 = lam3H3dUS + lam4H3dUS + lam5H3dUS;
      std::complex<double> Compile_338 = phi1*phi2*Compile_337;
      std::complex<double> Compile_339 = Compile_338 + Compile_51;
      std::complex<double> Compile_387 = -4.*Compile_12;
      std::complex<double> Compile_388 = 4.*Compile_13;
      std::complex<double> Compile_389 = -12.*lam1H3dUS*Compile_44;
      std::complex<double> Compile_390 = 2.*lam3H3dUS*Compile_44;
      std::complex<double> Compile_391 = 2.*lam4H3dUS*Compile_44;
      std::complex<double> Compile_392 = 2.*lam5H3dUS*Compile_44;
      std::complex<double> Compile_393 = 12.*lam2H3dUS*Compile_53;
      std::complex<double> Compile_394 = -2.*lam3H3dUS*Compile_53;
      std::complex<double> Compile_395 = -2.*lam4H3dUS*Compile_53;
      std::complex<double> Compile_396 = -2.*lam5H3dUS*Compile_53;
      std::complex<double> Compile_397 = pow(Compile_339,2);
      std::complex<double> Compile_398 = 64.*Compile_397;
      std::complex<double> Compile_399 = -6.*lam1H3dUS;
      std::complex<double> Compile_400 = 0. + lam3H3dUS + lam4H3dUS + \
      lam5H3dUS + Compile_399;
      std::complex<double> Compile_401 = 2.*Compile_400*Compile_44;
      std::complex<double> Compile_407 = 6.*lam2H3dUS;
      std::complex<double> Compile_408 = -1.*lam3H3dUS;
      std::complex<double> Compile_409 = -1.*lam4H3dUS;
      std::complex<double> Compile_410 = -1.*lam5H3dUS;
      std::complex<double> Compile_411 = Compile_407 + Compile_408 + \
      Compile_409 + Compile_410;
      std::complex<double> Compile_412 = 2.*Compile_411*Compile_53;
      std::complex<double> Compile_413 = 0. + Compile_387 + Compile_388 + \
      Compile_401 + Compile_412;
      std::complex<double> Compile_414 = pow(Compile_413,2);
      std::complex<double> Compile_415 = Compile_398 + Compile_414;
      std::complex<double> Compile_416 = sqrt(Compile_415);
      std::complex<double> Compile_417 = 0. + Compile_387 + Compile_388 + \
      Compile_389 + Compile_390 + Compile_391 + Compile_392 + Compile_393 + \
      Compile_394 + Compile_395 + Compile_396 + Compile_416;
      std::complex<double> Compile_340 = pow(Compile_339,-2);
      std::complex<double> Compile_418 = pow(Compile_417,2);
      std::complex<double> Compile_419 = 0.015625*Compile_340*Compile_418;
      std::complex<double> Compile_420 = 1. + Compile_419;
      std::complex<double> Compile_433 = pow(Compile_420,-1.5);
      std::complex<double> Compile_438 = 0. + lam3H3dUS + lam4H3dUS + \
      lam5H3dUS;
      std::complex<double> Compile_449 = 2.*lam1H3dUS;
      std::complex<double> Compile_450 = 0. + Compile_449;
      std::complex<double> Compile_439 = 1./Compile_339;
      std::complex<double> Compile_335 = 2.*lam2H3dUS;
      std::complex<double> Compile_336 = 0. + Compile_335;
      std::complex<double> Compile_452 = pow(Compile_339,-3);
      std::complex<double> Compile_453 = pow(Compile_417,3);
      std::complex<double> Compile_437 = -1.*phi2*Compile_336*Compile_433;
      std::complex<double> Compile_440 = \
      0.125*phi1*Compile_417*Compile_433*Compile_438*Compile_439;
      std::complex<double> Compile_441 = \
      -0.015625*phi2*Compile_340*Compile_418*Compile_433*Compile_438;
      std::complex<double> Compile_454 = \
      0.001953125*phi1*Compile_433*Compile_450*Compile_452*Compile_453;
      std::complex<double> Compile_455 = Compile_437 + Compile_440 + \
      Compile_441 + Compile_454;
      std::complex<double> Compile_456 = pow(Compile_455,2);
      std::complex<double> Compile_494 = 1./Compile_420;
      std::complex<double> Compile_458 = phi1*Compile_433*Compile_450;
      std::complex<double> Compile_459 = \
      0.125*phi2*Compile_417*Compile_433*Compile_438*Compile_439;
      std::complex<double> Compile_473 = \
      0.015625*phi1*Compile_340*Compile_418*Compile_433*Compile_438;
      std::complex<double> Compile_474 = \
      0.001953125*phi2*Compile_336*Compile_433*Compile_452*Compile_453;
      std::complex<double> Compile_475 = Compile_458 + Compile_459 + \
      Compile_473 + Compile_474;
      std::complex<double> Compile_482 = pow(Compile_475,2);
      std::complex<double> Compile_496 = 4.*Compile_12;
      std::complex<double> Compile_497 = 12.*lam1H3dUS*Compile_44;
      std::complex<double> Compile_521 = 2.*Compile_337*Compile_53;
      std::complex<double> Compile_522 = 0. + Compile_496 + Compile_497 + \
      Compile_521;
      std::complex<double> Compile_492 = 2.*Compile_438*Compile_44;
      std::complex<double> Compile_493 = 0. + Compile_388 + Compile_393 + \
      Compile_492;
      std::complex<double> Compile_524 = 0.125*Compile_417*Compile_439;
      std::complex<double> Compile_525 = tan(Compile_524);
      std::complex<double> Compile_526 = 2.*Compile_525;
      std::complex<double> Compile_527 = sin(Compile_526);
      std::complex<double> Compile_558 = 4.*Compile_51;
      std::complex<double> Compile_559 = lam4H3dUS + lam5H3dUS;
      std::complex<double> Compile_562 = 2.*phi1*phi2*Compile_559;
      std::complex<double> Compile_569 = Compile_558 + Compile_562;
      std::complex<double> Compile_556 = 4.*lam1H3dUS;
      std::complex<double> Compile_571 = -4.*Compile_13;
      std::complex<double> Compile_572 = 4.*lam1H3dUS*Compile_44;
      std::complex<double> Compile_573 = -2.*lam3H3dUS*Compile_44;
      std::complex<double> Compile_574 = -4.*lam2H3dUS*Compile_53;
      std::complex<double> Compile_575 = 2.*lam3H3dUS*Compile_53;
      std::complex<double> Compile_576 = pow(Compile_569,2);
      std::complex<double> Compile_580 = 4.*Compile_576;
      std::complex<double> Compile_584 = 4.*Compile_18;
      std::complex<double> Compile_585 = -2.*lam3H3dUS;
      std::complex<double> Compile_586 = 0. + Compile_556 + Compile_585;
      std::complex<double> Compile_587 = Compile_44*Compile_586;
      std::complex<double> Compile_590 = -4.*lam2H3dUS;
      std::complex<double> Compile_593 = 2.*lam3H3dUS;
      std::complex<double> Compile_594 = Compile_590 + Compile_593;
      std::complex<double> Compile_595 = Compile_53*Compile_594;
      std::complex<double> Compile_596 = 0. + Compile_584 + Compile_587 + \
      Compile_595;
      std::complex<double> Compile_597 = pow(Compile_596,2);
      std::complex<double> Compile_598 = Compile_580 + Compile_597;
      std::complex<double> Compile_599 = sqrt(Compile_598);
      std::complex<double> Compile_600 = 0. + Compile_496 + Compile_571 + \
      Compile_572 + Compile_573 + Compile_574 + Compile_575 + Compile_599;
      std::complex<double> Compile_570 = pow(Compile_569,-2);
      std::complex<double> Compile_601 = pow(Compile_600,2);
      std::complex<double> Compile_602 = 0.25*Compile_570*Compile_601;
      std::complex<double> Compile_604 = 1. + Compile_602;
      std::complex<double> Compile_608 = 1./Compile_604;
      std::complex<double> Compile_555 = 1./sqrt(Compile_420);
      std::complex<double> Compile_614 = 0. + Compile_593;
      std::complex<double> Compile_610 = 2.*Compile_559;
      std::complex<double> Compile_611 = 0. + Compile_610;
      std::complex<double> Compile_612 = 1./Compile_569;
      std::complex<double> Compile_557 = 0. + Compile_556;
      std::complex<double> Compile_609 = phi1*Compile_557*Compile_608;
      std::complex<double> Compile_613 = \
      -0.5*phi2*Compile_600*Compile_608*Compile_611*Compile_612;
      std::complex<double> Compile_615 = \
      0.25*phi1*Compile_570*Compile_601*Compile_608*Compile_614;
      std::complex<double> Compile_616 = Compile_609 + Compile_613 + \
      Compile_615;
      std::complex<double> Compile_618 = phi2*Compile_608*Compile_614;
      std::complex<double> Compile_619 = \
      -0.5*phi1*Compile_600*Compile_608*Compile_611*Compile_612;
      std::complex<double> Compile_620 = 4.*lam2H3dUS;
      std::complex<double> Compile_621 = 0. + Compile_620;
      std::complex<double> Compile_622 = \
      0.25*phi2*Compile_570*Compile_601*Compile_608*Compile_621;
      std::complex<double> Compile_624 = Compile_618 + Compile_619 + \
      Compile_622;
      std::complex<double> Compile_639 = phi1*Compile_608*Compile_614;
      std::complex<double> Compile_640 = \
      0.5*phi2*Compile_600*Compile_608*Compile_611*Compile_612;
      std::complex<double> Compile_641 = \
      0.25*phi1*Compile_557*Compile_570*Compile_601*Compile_608;
      std::complex<double> Compile_642 = Compile_639 + Compile_640 + \
      Compile_641;
      std::complex<double> Compile_644 = -1.*phi2*Compile_608*Compile_621;
      std::complex<double> Compile_647 = \
      -0.25*phi2*Compile_570*Compile_601*Compile_608*Compile_614;
      std::complex<double> Compile_648 = Compile_619 + Compile_644 + \
      Compile_647;
      std::complex<double> Compile_617 = \
      0.125*Compile_417*Compile_439*Compile_555*Compile_616;
      std::complex<double> Compile_626 = -1.*Compile_555*Compile_624;
      std::complex<double> Compile_627 = Compile_617 + Compile_626;
      std::complex<double> Compile_628 = pow(Compile_627,2);
      std::complex<double> Compile_495 = Compile_493*Compile_494;
      std::complex<double> Compile_523 = \
      0.015625*Compile_340*Compile_418*Compile_494*Compile_522;
      std::complex<double> Compile_533 = -4.*Compile_339*Compile_527;
      std::complex<double> Compile_534 = Compile_495 + Compile_523 + \
      Compile_533;
      std::complex<double> Compile_630 = -1.*Compile_555*Compile_616;
      std::complex<double> Compile_631 = \
      -0.125*Compile_417*Compile_439*Compile_555*Compile_624;
      std::complex<double> Compile_632 = Compile_630 + Compile_631;
      std::complex<double> Compile_633 = pow(Compile_632,2);
      std::complex<double> Compile_660 = Compile_44*Compile_557;
      std::complex<double> Compile_666 = 0. + Compile_496 + Compile_575 + \
      Compile_660;
      std::complex<double> Compile_667 = Compile_608*Compile_666;
      std::complex<double> Compile_668 = Compile_44*Compile_614;
      std::complex<double> Compile_669 = 4.*lam2H3dUS*Compile_53;
      std::complex<double> Compile_670 = 0. + Compile_388 + Compile_668 + \
      Compile_669;
      std::complex<double> Compile_671 = \
      0.25*Compile_570*Compile_601*Compile_608*Compile_670;
      std::complex<double> Compile_672 = 0.5*Compile_600*Compile_612;
      std::complex<double> Compile_673 = tan(Compile_672);
      std::complex<double> Compile_674 = 2.*Compile_673;
      std::complex<double> Compile_683 = sin(Compile_674);
      std::complex<double> Compile_688 = -1.*Compile_569*Compile_683;
      std::complex<double> Compile_692 = Compile_667 + Compile_671 + \
      Compile_688;
      std::complex<double> Compile_693 = sqrt(Compile_692);
      std::complex<double> Compile_694 = sqrt(Compile_534);
      std::complex<double> Compile_695 = 0.5*Compile_694;
      std::complex<double> Compile_696 = Compile_693 + Compile_695;
      std::complex<double> Compile_698 = 1./Compile_696;
      std::complex<double> Compile_699 = mu3US*Compile_698;
      std::complex<double> Compile_700 = log(Compile_699);
      std::complex<double> Compile_643 = \
      0.125*Compile_417*Compile_439*Compile_555*Compile_642;
      std::complex<double> Compile_649 = Compile_555*Compile_648;
      std::complex<double> Compile_650 = Compile_643 + Compile_649;
      std::complex<double> Compile_651 = pow(Compile_650,2);
      std::complex<double> Compile_653 = -1.*Compile_555*Compile_642;
      std::complex<double> Compile_655 = \
      0.125*Compile_417*Compile_439*Compile_555*Compile_648;
      std::complex<double> Compile_656 = Compile_653 + Compile_655;
      std::complex<double> Compile_657 = pow(Compile_656,2);
      std::complex<double> Compile_720 = Compile_608*Compile_670;
      std::complex<double> Compile_721 = \
      0.25*Compile_570*Compile_601*Compile_608*Compile_666;
      std::complex<double> Compile_722 = Compile_569*Compile_683;
      std::complex<double> Compile_723 = Compile_720 + Compile_721 + \
      Compile_722;
      std::complex<double> Compile_724 = sqrt(Compile_723);
      std::complex<double> Compile_543 = Compile_494*Compile_522;
      std::complex<double> Compile_544 = \
      0.015625*Compile_340*Compile_418*Compile_493*Compile_494;
      std::complex<double> Compile_545 = 4.*Compile_339*Compile_527;
      std::complex<double> Compile_546 = Compile_543 + Compile_544 + \
      Compile_545;
      std::complex<double> Compile_744 = 2.*Compile_51;
      std::complex<double> Compile_745 = 2.*lam5H3dUS*phi1*phi2;
      std::complex<double> Compile_746 = Compile_744 + Compile_745;
      std::complex<double> Compile_748 = -2.*lam4H3dUS*Compile_44;
      std::complex<double> Compile_750 = 2.*lam4H3dUS*Compile_53;
      std::complex<double> Compile_751 = pow(Compile_746,2);
      std::complex<double> Compile_752 = 16.*Compile_751;
      std::complex<double> Compile_753 = 0. + lam5H3dUS + Compile_408 + \
      Compile_409 + Compile_449;
      std::complex<double> Compile_754 = -2.*Compile_44*Compile_753;
      std::complex<double> Compile_755 = lam5H3dUS + Compile_335 + \
      Compile_408 + Compile_409;
      std::complex<double> Compile_756 = 2.*Compile_53*Compile_755;
      std::complex<double> Compile_757 = 0. + Compile_387 + Compile_388 + \
      Compile_754 + Compile_756;
      std::complex<double> Compile_758 = pow(Compile_757,2);
      std::complex<double> Compile_759 = Compile_752 + Compile_758;
      std::complex<double> Compile_760 = sqrt(Compile_759);
      std::complex<double> Compile_761 = 0. + Compile_392 + Compile_396 + \
      Compile_496 + Compile_571 + Compile_572 + Compile_573 + Compile_574 + \
      Compile_575 + Compile_748 + Compile_750 + Compile_760;
      std::complex<double> Compile_747 = pow(Compile_746,-2);
      std::complex<double> Compile_762 = pow(Compile_761,2);
      std::complex<double> Compile_763 = 0.0625*Compile_747*Compile_762;
      std::complex<double> Compile_764 = 1. + Compile_763;
      std::complex<double> Compile_765 = 1./Compile_764;
      std::complex<double> Compile_771 = 0. + lam3H3dUS + lam4H3dUS + \
      Compile_410;
      std::complex<double> Compile_769 = 1./Compile_746;
      std::complex<double> Compile_776 = 1./sqrt(Compile_764);
      std::complex<double> Compile_792 = lam3H3dUS + lam4H3dUS + Compile_410;
      std::complex<double> Compile_767 = 2.*lam5H3dUS;
      std::complex<double> Compile_768 = 0. + Compile_767;
      std::complex<double> Compile_876 = phi1*Compile_765*Compile_771;
      std::complex<double> Compile_877 = \
      0.25*phi2*Compile_761*Compile_765*Compile_768*Compile_769;
      std::complex<double> Compile_878 = \
      0.0625*phi1*Compile_450*Compile_747*Compile_762*Compile_765;
      std::complex<double> Compile_879 = Compile_876 + Compile_877 + \
      Compile_878;
      std::complex<double> Compile_894 = -1.*phi2*Compile_336*Compile_765;
      std::complex<double> Compile_895 = \
      -0.25*phi1*Compile_761*Compile_765*Compile_768*Compile_769;
      std::complex<double> Compile_896 = \
      -0.0625*phi2*Compile_747*Compile_762*Compile_765*Compile_771;
      std::complex<double> Compile_897 = Compile_894 + Compile_895 + \
      Compile_896;
      std::complex<double> Compile_766 = phi1*Compile_450*Compile_765;
      std::complex<double> Compile_770 = \
      -0.25*phi2*Compile_761*Compile_765*Compile_768*Compile_769;
      std::complex<double> Compile_772 = \
      0.0625*phi1*Compile_747*Compile_762*Compile_765*Compile_771;
      std::complex<double> Compile_773 = Compile_766 + Compile_770 + \
      Compile_772;
      std::complex<double> Compile_775 = phi2*Compile_765*Compile_771;
      std::complex<double> Compile_777 = 2.*lam5H3dUS*phi1*Compile_776;
      std::complex<double> Compile_778 = \
      -0.25*phi2*Compile_336*Compile_761*Compile_769*Compile_776;
      std::complex<double> Compile_779 = Compile_777 + Compile_778;
      std::complex<double> Compile_785 = \
      -0.25*Compile_761*Compile_769*Compile_776*Compile_779;
      std::complex<double> Compile_786 = 0. + Compile_775 + Compile_785;
      std::complex<double> Compile_844 = -1.*phi1*Compile_433*Compile_438;
      std::complex<double> Compile_846 = -3.*lam2H3dUS;
      std::complex<double> Compile_847 = lam3H3dUS + lam4H3dUS + lam5H3dUS \
      + Compile_846;
      std::complex<double> Compile_848 = 2.*Compile_847;
      std::complex<double> Compile_849 = 0. + Compile_848;
      std::complex<double> Compile_856 = \
      0.125*phi2*Compile_417*Compile_433*Compile_439*Compile_849;
      std::complex<double> Compile_857 = \
      -0.046875*phi1*Compile_340*Compile_418*Compile_433*Compile_450;
      std::complex<double> Compile_858 = -2.*phi1*Compile_337*Compile_555;
      std::complex<double> Compile_859 = \
      0.125*phi2*Compile_417*Compile_438*Compile_439*Compile_555;
      std::complex<double> Compile_860 = Compile_858 + Compile_859;
      std::complex<double> Compile_861 = \
      -0.015625*Compile_340*Compile_418*Compile_494*Compile_860;
      std::complex<double> Compile_874 = Compile_844 + Compile_856 + \
      Compile_857 + Compile_861;
      std::complex<double> Compile_875 = pow(Compile_874,2);
      std::complex<double> Compile_734 = sqrt(Compile_546);
      std::complex<double> Compile_735 = 0.5*Compile_734;
      std::complex<double> Compile_825 = -1.*phi2*Compile_433*Compile_438;
      std::complex<double> Compile_826 = -3.*lam1H3dUS;
      std::complex<double> Compile_827 = lam3H3dUS + lam4H3dUS + lam5H3dUS \
      + Compile_826;
      std::complex<double> Compile_828 = 2.*Compile_827;
      std::complex<double> Compile_829 = 0. + Compile_828;
      std::complex<double> Compile_835 = \
      -0.125*phi1*Compile_417*Compile_433*Compile_439*Compile_829;
      std::complex<double> Compile_836 = \
      -0.046875*phi2*Compile_336*Compile_340*Compile_418*Compile_433;
      std::complex<double> Compile_837 = 2.*phi2*Compile_337*Compile_555;
      std::complex<double> Compile_838 = \
      0.125*phi1*Compile_417*Compile_438*Compile_439*Compile_555;
      std::complex<double> Compile_839 = Compile_837 + Compile_838;
      std::complex<double> Compile_841 = \
      0.015625*Compile_340*Compile_418*Compile_494*Compile_839;
      std::complex<double> Compile_842 = Compile_825 + Compile_835 + \
      Compile_836 + Compile_841;
      std::complex<double> Compile_843 = pow(Compile_842,2);
      std::complex<double> Compile_908 = -1.*Compile_555*Compile_773;
      std::complex<double> Compile_909 = \
      -0.125*Compile_417*Compile_439*Compile_555*Compile_786;
      std::complex<double> Compile_910 = Compile_908 + Compile_909;
      std::complex<double> Compile_922 = pow(Compile_910,2);
      std::complex<double> Compile_793 = 2.*Compile_53*Compile_792;
      std::complex<double> Compile_794 = 0. + Compile_496 + Compile_572 + \
      Compile_793;
      std::complex<double> Compile_795 = Compile_765*Compile_794;
      std::complex<double> Compile_796 = 2.*Compile_44*Compile_792;
      std::complex<double> Compile_797 = 0. + Compile_388 + Compile_669 + \
      Compile_796;
      std::complex<double> Compile_798 = \
      0.0625*Compile_747*Compile_762*Compile_765*Compile_797;
      std::complex<double> Compile_799 = 0.25*Compile_761*Compile_769;
      std::complex<double> Compile_800 = tan(Compile_799);
      std::complex<double> Compile_806 = 2.*Compile_800;
      std::complex<double> Compile_807 = sin(Compile_806);
      std::complex<double> Compile_808 = -2.*Compile_746*Compile_807;
      std::complex<double> Compile_815 = Compile_795 + Compile_798 + \
      Compile_808;
      std::complex<double> Compile_816 = sqrt(Compile_815);
      std::complex<double> Compile_882 = \
      0.125*Compile_417*Compile_439*Compile_555*Compile_879;
      std::complex<double> Compile_898 = Compile_555*Compile_897;
      std::complex<double> Compile_902 = Compile_882 + Compile_898;
      std::complex<double> Compile_903 = pow(Compile_902,2);
      std::complex<double> Compile_904 = -1.*Compile_555*Compile_879;
      std::complex<double> Compile_905 = \
      0.125*Compile_417*Compile_439*Compile_555*Compile_897;
      std::complex<double> Compile_906 = Compile_904 + Compile_905;
      std::complex<double> Compile_907 = pow(Compile_906,2);
      std::complex<double> Compile_940 = Compile_765*Compile_797;
      std::complex<double> Compile_941 = \
      0.0625*Compile_747*Compile_762*Compile_765*Compile_794;
      std::complex<double> Compile_942 = 2.*Compile_746*Compile_807;
      std::complex<double> Compile_948 = Compile_940 + Compile_941 + \
      Compile_942;
      std::complex<double> Compile_949 = sqrt(Compile_948);
      std::complex<double> Compile_31 = Compile_12 + Compile_13 + Compile_24;
      std::complex<double> Compile_969 = Compile_12 + Compile_17 + \
      Compile_24;
      std::complex<double> Compile_971 = pow(m12R3dUS,-4);
      std::complex<double> Compile_972 = pow(Compile_969,2);
      std::complex<double> Compile_973 = 0.25*Compile_971*Compile_972;
      std::complex<double> Compile_974 = 1. + Compile_973;
      std::complex<double> Compile_975 = pow(Compile_974,-2);
      std::complex<double> Compile_26 = -1.*Compile_24;
      std::complex<double> Compile_27 = Compile_12 + Compile_13 + Compile_26;
      std::complex<double> Compile_968 = pow(m12R3dUS,-8);
      std::complex<double> Compile_970 = pow(Compile_969,4);
      std::complex<double> Compile_978 = \
      0.5*Compile_337*Compile_971*Compile_972*Compile_975;
      std::complex<double> Compile_976 = \
      0.125*lam1H3dUS*Compile_968*Compile_970*Compile_975;
      std::complex<double> Compile_977 = 2.*lam2H3dUS*Compile_975;
      std::complex<double> Compile_979 = Compile_976 + Compile_977 + \
      Compile_978;
      std::complex<double> Compile_981 = 2.*lam1H3dUS*Compile_975;
      std::complex<double> Compile_982 = \
      0.125*lam2H3dUS*Compile_968*Compile_970*Compile_975;
      std::complex<double> Compile_983 = Compile_978 + Compile_981 + \
      Compile_982;
      std::complex<double> Compile_989 = -1.*Compile_12;
      std::complex<double> Compile_990 = Compile_13 + Compile_24 + \
      Compile_989;
      std::complex<double> Compile_992 = pow(Compile_990,2);
      std::complex<double> Compile_993 = 0.25*Compile_971*Compile_992;
      std::complex<double> Compile_994 = 1. + Compile_993;
      std::complex<double> Compile_995 = pow(Compile_994,-2);
      std::complex<double> Compile_991 = pow(Compile_990,4);
      std::complex<double> Compile_998 = \
      0.5*Compile_337*Compile_971*Compile_992*Compile_995;
      std::complex<double> Compile_1027 = 1./sqrt(Compile_27);
      std::complex<double> Compile_1028 = \
      0.7071067811865475*mu3US*Compile_1027;
      std::complex<double> Compile_1029 = log(Compile_1028);
      std::complex<double> Compile_1037 = 1./sqrt(Compile_31);
      std::complex<double> Compile_1038 = \
      0.7071067811865475*mu3US*Compile_1037;
      std::complex<double> Compile_1039 = log(Compile_1038);
      std::complex<double> Compile_1016 = pow(m12R3dUS,-2);
      std::complex<double> Compile_1017 = 0.5*Compile_1016*Compile_969;
      std::complex<double> Compile_1018 = tan(Compile_1017);
      std::complex<double> Compile_1008 = sqrt(Compile_27);
      std::complex<double> Compile_1009 = sqrt(Compile_31);
      std::complex<double> Compile_1049 = 0.5*Compile_1016*Compile_990;
      std::complex<double> Compile_1050 = tan(Compile_1049);
      std::complex<double> Compile_1051 = Compile_1018 + Compile_1050;
      std::complex<double> Compile_1052 = cos(Compile_1051);
      std::complex<double> Compile_1053 = pow(Compile_1052,2);
      std::complex<double> Compile_1061 = 0.5*Compile_27;
      std::complex<double> Compile_1062 = 0.5*Compile_31;
      std::complex<double> Compile_1063 = Compile_1061 + Compile_1062;
      std::complex<double> Compile_1065 = 0.7071067811865475*Compile_1008;
      std::complex<double> Compile_1069 = 0.7071067811865475*Compile_1009;
      std::complex<double> Compile_1070 = Compile_1065 + Compile_1069;
      std::complex<double> Compile_1071 = 1./Compile_1070;
      std::complex<double> Compile_1072 = mu3US*Compile_1071;
      std::complex<double> Compile_1073 = log(Compile_1072);
      std::complex<double> Compile_1088 = 1./Compile_974;
      std::complex<double> Compile_1087 = 1./Compile_994;
      std::complex<double> Compile_1090 = Compile_1088*Compile_792;
      std::complex<double> Compile_1101 = \
      0.25*Compile_1088*Compile_792*Compile_971*Compile_972;
      std::complex<double> Compile_1104 = \
      lam5H3dUS*Compile_1087*Compile_1088*Compile_969*Compile_971*Compile_990;
      std::complex<double> Compile_1012 = -2.*lam4H3dUS;
      std::complex<double> Compile_1013 = -1.*lam1H3dUS;
      std::complex<double> Compile_1014 = -1.*lam2H3dUS;
      std::complex<double> Compile_1015 = lam3H3dUS + lam4H3dUS + lam5H3dUS \
      + Compile_1013 + Compile_1014;
      std::complex<double> Compile_1019 = 4.*Compile_1018;
      std::complex<double> Compile_1020 = cos(Compile_1019);
      std::complex<double> Compile_1119 = -6.*lam2H3dUS;
      std::complex<double> Compile_1121 = -2.*lam5H3dUS;
      std::complex<double> Compile_1110 = \
      0.5*lam1H3dUS*Compile_1088*Compile_971*Compile_972;
      std::complex<double> Compile_1112 = Compile_1090 + Compile_1110;
      std::complex<double> Compile_1135 = \
      -0.0625*lam3H3dUS*Compile_1087*Compile_1088*Compile_968*Compile_972*\
      Compile_992;
      std::complex<double> Compile_1136 = \
      -0.0625*lam4H3dUS*Compile_1087*Compile_1088*Compile_968*Compile_972*\
      Compile_992;
      std::complex<double> Compile_1089 = \
      0.5*lam2H3dUS*Compile_1088*Compile_971*Compile_972;
      std::complex<double> Compile_1091 = Compile_1089 + Compile_1090;
      std::complex<double> Compile_1138 = \
      0.0625*lam5H3dUS*Compile_1087*Compile_1088*Compile_968*Compile_972*\
      Compile_992;
      std::complex<double> Compile_1139 = 2.*Compile_1018;
      std::complex<double> Compile_1140 = sin(Compile_1139);
      std::complex<double> Compile_1141 = 2.*Compile_1050;
      std::complex<double> Compile_1146 = sin(Compile_1141);
      std::complex<double> Compile_1147 = \
      -1.*lam5H3dUS*Compile_1140*Compile_1146;
      std::complex<double> Compile_1157 = 2.*lam3H3dUS*Compile_1088;
      std::complex<double> Compile_1164 = \
      0.5*lam3H3dUS*Compile_1088*Compile_971*Compile_972;
      std::complex<double> Compile_1156 = \
      lam1H3dUS*Compile_1088*Compile_971*Compile_972;
      std::complex<double> Compile_1159 = Compile_1156 + Compile_1157;
      std::complex<double> Compile_1163 = 4.*lam2H3dUS*Compile_1088;
      std::complex<double> Compile_1165 = Compile_1163 + Compile_1164;
      std::complex<double> Compile_1193 = pow(Compile_1140,2);
      std::complex<double> Compile_1189 = \
      -0.125*lam3H3dUS*Compile_968*Compile_970*Compile_975;
      std::complex<double> Compile_1170 = \
      lam2H3dUS*Compile_1088*Compile_971*Compile_972;
      std::complex<double> Compile_1171 = Compile_1157 + Compile_1170;
      std::complex<double> Compile_1194 = lam4H3dUS*Compile_1193;
      std::complex<double> Compile_1195 = lam5H3dUS*Compile_1193;
      std::complex<double> Compile_1173 = 4.*lam1H3dUS*Compile_1088;
      std::complex<double> Compile_1174 = Compile_1164 + Compile_1173;
      std::complex<double> Compile_1215 = \
      -0.125*lam3H3dUS*Compile_1087*Compile_1088*Compile_968*Compile_972*\
      Compile_992;
      std::complex<double> Compile_1217 = \
      -1.*lam4H3dUS*Compile_1140*Compile_1146;
      std::complex<double> Compile_1024 = 1.5*Compile_27;
      std::complex<double> Compile_1031 = 2.*Compile_1029*Compile_27;
      std::complex<double> Compile_1032 = Compile_1024 + Compile_1031;
      std::complex<double> Compile_1034 = -2.*Compile_31;
      std::complex<double> Compile_1035 = Compile_1034 + Compile_17 + \
      Compile_26 + Compile_989;
      std::complex<double> Compile_1036 = 0.5*Compile_1035;
      std::complex<double> Compile_1040 = -2.*Compile_1039*Compile_31;
      std::complex<double> Compile_1041 = Compile_1036 + Compile_1040;
      std::complex<double> Compile_1243 = sin(Compile_1051);
      std::complex<double> Compile_1244 = pow(Compile_1243,2);
      std::complex<double> Compile_1054 = Compile_17 + Compile_26 + \
      Compile_989;
      std::complex<double> Compile_1055 = 0.5*Compile_1054;
      std::complex<double> Compile_1056 = Compile_17 + Compile_24 + \
      Compile_989;
      std::complex<double> Compile_1057 = 0.5*Compile_1056;
      std::complex<double> Compile_1058 = -2.*Compile_1008*Compile_1009;
      std::complex<double> Compile_1059 = Compile_1055 + Compile_1057 + \
      Compile_1058;
      std::complex<double> Compile_1060 = 0.5*Compile_1059;
      std::complex<double> Compile_1074 = -2.*Compile_1063*Compile_1073;
      std::complex<double> Compile_1075 = Compile_1060 + Compile_1074;
      std::complex<double> Compile_1245 = \
      0.0015831434944115277*Compile_1032*Compile_1244*Compile_64;
      std::complex<double> Compile_1246 = \
      -0.0015831434944115277*Compile_1041*Compile_1244*Compile_64;
      std::complex<double> Compile_1010 = -2.*lam1H3dUS;
      std::complex<double> Compile_1011 = -2.*lam2H3dUS;
      std::complex<double> Compile_65 = Compile_61*Compile_64;
      std::complex<double> Compile_1270 = sqrt(Compile_65);
      std::complex<double> Compile_1279 = pow(Compile_604,-2);
      std::complex<double> Compile_1281 = 2.*Compile_337;
      std::complex<double> Compile_1282 = 0. + Compile_1281;
      std::complex<double> Compile_1283 = \
      0.25*Compile_1279*Compile_1282*Compile_570*Compile_601;
      std::complex<double> Compile_1284 = pow(Compile_569,-4);
      std::complex<double> Compile_1288 = pow(Compile_600,4);
      std::complex<double> Compile_1298 = lam4H3dUS + Compile_410;
      std::complex<double> Compile_1299 = 2.*Compile_1298;
      std::complex<double> Compile_1300 = 0. + Compile_1299;
      std::complex<double> Compile_1301 = pow(Compile_1300,2);
      std::complex<double> Compile_1302 = -1.*phi2*Compile_776;
      std::complex<double> Compile_1303 = \
      -0.25*phi1*Compile_761*Compile_769*Compile_776;
      std::complex<double> Compile_1304 = Compile_1302 + Compile_1303;
      std::complex<double> Compile_1305 = pow(Compile_1304,2);
      std::complex<double> Compile_1307 = phi1*Compile_776;
      std::complex<double> Compile_1308 = \
      -0.25*phi2*Compile_761*Compile_769*Compile_776;
      std::complex<double> Compile_1309 = Compile_1307 + Compile_1308;
      std::complex<double> Compile_1310 = pow(Compile_1309,2);
      std::complex<double> Compile_1312 = 0.5*Compile_693;
      std::complex<double> Compile_1317 = 0.5*Compile_724;
      std::complex<double> Compile_1331 = cos(Compile_674);
      std::complex<double> Compile_1332 = phi2*Compile_1331*Compile_611;
      std::complex<double> Compile_1333 = 0. + lam3H3dUS + Compile_1010;
      std::complex<double> Compile_1335 = 2.*Compile_1333;
      std::complex<double> Compile_1336 = 0. + Compile_1335;
      std::complex<double> Compile_1337 = \
      -1.*phi1*Compile_1336*Compile_683;
      std::complex<double> Compile_1338 = Compile_1332 + Compile_1337;
      std::complex<double> Compile_1339 = \
      0.125*Compile_1338*Compile_417*Compile_439*Compile_555;
      std::complex<double> Compile_1340 = phi1*Compile_1331*Compile_611;
      std::complex<double> Compile_1341 = 0. + lam3H3dUS;
      std::complex<double> Compile_1342 = -2.*Compile_1341;
      std::complex<double> Compile_1343 = 0. + Compile_1342 + Compile_620;
      std::complex<double> Compile_1344 = \
      -1.*phi2*Compile_1343*Compile_683;
      std::complex<double> Compile_1345 = Compile_1340 + Compile_1344;
      std::complex<double> Compile_1347 = -1.*Compile_1345*Compile_555;
      std::complex<double> Compile_1348 = Compile_1339 + Compile_1347;
      std::complex<double> Compile_1349 = pow(Compile_1348,2);
      std::complex<double> Compile_1356 = -1.*Compile_1338*Compile_555;
      std::complex<double> Compile_1357 = \
      -0.125*Compile_1345*Compile_417*Compile_439*Compile_555;
      std::complex<double> Compile_1358 = Compile_1356 + Compile_1357;
      std::complex<double> Compile_1359 = pow(Compile_1358,2);
      std::complex<double> Compile_273 = sqrt(Compile_64);
      std::complex<double> Compile_274 = 0.5*Compile_269*Compile_273;
      std::complex<double> Compile_1385 = Compile_269*Compile_273;
      std::complex<double> Compile_1368 = 1./Compile_61;
      std::complex<double> Compile_1422 = 0.5*g23dUS*Compile_269;
      std::complex<double> Compile_1431 = -0.5*Compile_724;
      std::complex<double> Compile_1429 = -0.5*g23dUS*Compile_269;
      std::complex<double> Compile_1424 = Compile_1317 + Compile_1422 + \
      Compile_695;
      std::complex<double> Compile_1450 = 0.25*Compile_61*Compile_64;
      std::complex<double> Compile_1442 = -1.*Compile_608*Compile_670;
      std::complex<double> Compile_1443 = \
      -0.25*Compile_570*Compile_601*Compile_608*Compile_666;
      std::complex<double> Compile_1444 = Compile_1442 + Compile_1443 + \
      Compile_688;
      std::complex<double> Compile_1445 = 0.25*Compile_1444;
      std::complex<double> Compile_1446 = 0.25*Compile_534;
      std::complex<double> Compile_1419 = Compile_525 + Compile_673;
      std::complex<double> Compile_1420 = cos(Compile_1419);
      std::complex<double> Compile_1421 = pow(Compile_1420,2);
      std::complex<double> Compile_1474 = 0.25*Compile_692;
      std::complex<double> Compile_1475 = -1.*Compile_494*Compile_522;
      std::complex<double> Compile_1476 = \
      -0.015625*Compile_340*Compile_418*Compile_493*Compile_494;
      std::complex<double> Compile_1477 = Compile_1475 + Compile_1476 + \
      Compile_533;
      std::complex<double> Compile_1478 = 0.25*Compile_1477;
      std::complex<double> Compile_1489 = -0.5*Compile_693;
      std::complex<double> Compile_1484 = Compile_1312 + Compile_1422 + \
      Compile_735;
      std::complex<double> Compile_547 = 1./sqrt(Compile_546);
      std::complex<double> Compile_49 = pow(phi1,4);
      std::complex<double> Compile_58 = pow(phi2,4);
      std::complex<double> Compile_270 = g23dUS*Compile_269;
      std::complex<double> Compile_1516 = pow(Compile_546,2);
      std::complex<double> Compile_1506 = pow(Compile_61,-2);
      std::complex<double> Compile_536 = 1./sqrt(Compile_534);
      std::complex<double> Compile_1452 = -1.*Compile_493*Compile_494;
      std::complex<double> Compile_1453 = \
      -0.015625*Compile_340*Compile_418*Compile_494*Compile_522;
      std::complex<double> Compile_1454 = Compile_1452 + Compile_1453 + \
      Compile_545;
      std::complex<double> Compile_1457 = 0.25*Compile_1454;
      std::complex<double> Compile_1533 = -0.5*Compile_49*Compile_60;
      std::complex<double> Compile_1534 = \
      -1.*Compile_44*Compile_53*Compile_60;
      std::complex<double> Compile_1535 = -0.5*Compile_58*Compile_60;
      std::complex<double> Compile_1536 = Compile_1533 + Compile_1534 + \
      Compile_1535;
      std::complex<double> Compile_1537 = 0.25*Compile_1536;
      std::complex<double> Compile_1543 = 0.5*Compile_49*Compile_60;
      std::complex<double> Compile_1544 = Compile_44*Compile_53*Compile_60;
      std::complex<double> Compile_1545 = 0.5*Compile_58*Compile_60;
      std::complex<double> Compile_1564 = pow(Compile_534,2);
      std::complex<double> Compile_1528 = 0.5*Compile_61*Compile_64;
      std::complex<double> Compile_1507 = phi1*Compile_555;
      std::complex<double> Compile_1508 = \
      0.125*phi2*Compile_417*Compile_439*Compile_555;
      std::complex<double> Compile_1509 = Compile_1507 + Compile_1508;
      std::complex<double> Compile_1510 = pow(Compile_1509,2);
      std::complex<double> Compile_1511 = \
      -0.00039578587360288194*Compile_1270*Compile_61*Compile_64*Compile_734;
      std::complex<double> Compile_1512 = 2.*mu3US*Compile_547;
      std::complex<double> Compile_1513 = log(Compile_1512);
      std::complex<double> Compile_1514 = 2.*Compile_1513;
      std::complex<double> Compile_1515 = 1. + Compile_1514;
      std::complex<double> Compile_1517 = \
      0.00019789293680144097*Compile_1515*Compile_1516;
      std::complex<double> Compile_1524 = Compile_1450 + Compile_1478;
      std::complex<double> Compile_1525 = pow(Compile_1524,2);
      std::complex<double> Compile_1527 = \
      0.009947183943243459*Compile_61*Compile_64*Compile_734;
      std::complex<double> Compile_1529 = Compile_1478 + Compile_1528;
      std::complex<double> Compile_1530 = \
      -0.039788735772973836*Compile_1270*Compile_1529;
      std::complex<double> Compile_1531 = Compile_1527 + Compile_1530;
      std::complex<double> Compile_1532 = \
      -0.039788735772973836*Compile_1270*Compile_1531;
      std::complex<double> Compile_1603 = pow(Compile_61,2);
      std::complex<double> Compile_1555 = phi2*Compile_555;
      std::complex<double> Compile_1556 = \
      -0.125*phi1*Compile_417*Compile_439*Compile_555;
      std::complex<double> Compile_1557 = Compile_1555 + Compile_1556;
      std::complex<double> Compile_1558 = pow(Compile_1557,2);
      std::complex<double> Compile_1559 = \
      -0.00039578587360288194*Compile_1270*Compile_61*Compile_64*Compile_694;
      std::complex<double> Compile_1560 = 2.*mu3US*Compile_536;
      std::complex<double> Compile_1561 = log(Compile_1560);
      std::complex<double> Compile_1562 = 2.*Compile_1561;
      std::complex<double> Compile_1563 = 1. + Compile_1562;
      std::complex<double> Compile_1565 = \
      0.00019789293680144097*Compile_1563*Compile_1564;
      std::complex<double> Compile_1572 = Compile_1450 + Compile_1457;
      std::complex<double> Compile_1573 = pow(Compile_1572,2);
      std::complex<double> Compile_1587 = \
      0.009947183943243459*Compile_61*Compile_64*Compile_694;
      std::complex<double> Compile_1588 = Compile_1457 + Compile_1528;
      std::complex<double> Compile_1589 = \
      -0.039788735772973836*Compile_1270*Compile_1588;
      std::complex<double> Compile_1590 = Compile_1587 + Compile_1589;
      std::complex<double> Compile_1591 = \
      -0.039788735772973836*Compile_1270*Compile_1590;
      std::complex<double> Compile_1604 = -0.125*Compile_1603*Compile_60;
      std::complex<double> Compile_1610 = 0.0625*Compile_1603*Compile_60;
      std::complex<double> Compile_1611 = 1.5*Compile_61*Compile_64;
      std::complex<double> Compile_1644 = Compile_1312 + Compile_1422 + \
      Compile_695;
      std::complex<double> Compile_1494 = \
      0.009947183943243459*Compile_61*Compile_64*Compile_693;
      std::complex<double> Compile_1495 = -1.*Compile_608*Compile_666;
      std::complex<double> Compile_1496 = \
      -0.25*Compile_570*Compile_601*Compile_608*Compile_670;
      std::complex<double> Compile_1497 = Compile_1495 + Compile_1496 + \
      Compile_722;
      std::complex<double> Compile_1498 = 0.25*Compile_1497;
      std::complex<double> Compile_1451 = 0.25*Compile_723;
      std::complex<double> Compile_1674 = Compile_1317 + Compile_1422 + \
      Compile_735;
      std::complex<double> Compile_1499 = 0.25*Compile_546;
      std::complex<double> Compile_1461 = \
      0.009947183943243459*Compile_61*Compile_64*Compile_724;
      std::complex<double> Compile_1669 = sin(Compile_1419);
      std::complex<double> Compile_1670 = pow(Compile_1669,2);
      std::complex<double> Compile_1708 = pow(Compile_420,-2);
      std::complex<double> Compile_1710 = \
      0.03125*Compile_1708*Compile_340*Compile_418*Compile_438;
      std::complex<double> Compile_1711 = pow(Compile_339,-4);
      std::complex<double> Compile_1712 = pow(Compile_417,4);
      std::complex<double> Compile_1720 = pow(Compile_764,-2);
      std::complex<double> Compile_1722 = \
      0.125*Compile_1720*Compile_337*Compile_747*Compile_762;
      std::complex<double> Compile_1723 = pow(Compile_746,-4);
      std::complex<double> Compile_1724 = pow(Compile_761,4);
      std::complex<double> Compile_1737 = 2.*Compile_1341;
      std::complex<double> Compile_1738 = 0. + Compile_1737;
      std::complex<double> Compile_1739 = Compile_1738*Compile_608;
      std::complex<double> Compile_1744 = \
      0.25*Compile_1738*Compile_570*Compile_601*Compile_608;
      std::complex<double> Compile_1755 = 1./sqrt(Compile_604);
      std::complex<double> Compile_1740 = \
      0.25*Compile_557*Compile_570*Compile_601*Compile_608;
      std::complex<double> Compile_1741 = 0. + Compile_1739 + Compile_1740;
      std::complex<double> Compile_1764 = \
      -0.0078125*lam3H3dUS*Compile_340*Compile_418*Compile_494*Compile_570*\
      Compile_601*Compile_608;
      std::complex<double> Compile_1749 = \
      0.25*Compile_570*Compile_601*Compile_608*Compile_621;
      std::complex<double> Compile_1750 = 0. + Compile_1739 + Compile_1749;
      std::complex<double> Compile_1766 = \
      -1.*lam4H3dUS*Compile_527*Compile_683;
      std::complex<double> Compile_1767 = \
      -1.*lam5H3dUS*Compile_527*Compile_683;
      std::complex<double> Compile_1752 = Compile_557*Compile_608;
      std::complex<double> Compile_1774 = Compile_608*Compile_614;
      std::complex<double> Compile_1743 = Compile_608*Compile_621;
      std::complex<double> Compile_1775 = 0. + Compile_1749 + Compile_1774;
      std::complex<double> Compile_1802 = \
      -0.03125*lam3H3dUS*Compile_570*Compile_601*Compile_608*Compile_747*\
      Compile_762*Compile_765;
      std::complex<double> Compile_1786 = \
      -0.5*Compile_1755*Compile_557*Compile_600*Compile_612;
      std::complex<double> Compile_1787 = 0. + Compile_1786;
      std::complex<double> Compile_1788 = \
      -0.5*Compile_1755*Compile_1787*Compile_600*Compile_612;
      std::complex<double> Compile_1789 = Compile_1774 + Compile_1788;
      std::complex<double> Compile_1804 = \
      lam4H3dUS*Compile_683*Compile_807;
      std::complex<double> Compile_1805 = \
      lam5H3dUS*Compile_683*Compile_807;
      std::complex<double> Compile_1260 = 2.*Compile_1015;
      std::complex<double> Compile_1261 = 0. + Compile_1260;
      std::complex<double> Compile_1822 = Compile_765*Compile_771;
      std::complex<double> Compile_1827 = \
      0.0625*Compile_747*Compile_762*Compile_765*Compile_771;
      std::complex<double> Compile_1835 = \
      0.0625*Compile_450*Compile_747*Compile_762*Compile_765;
      std::complex<double> Compile_1836 = 0. + Compile_1822 + Compile_1835;
      std::complex<double> Compile_1849 = \
      -0.0009765625*lam3H3dUS*Compile_340*Compile_418*Compile_494*Compile_747*Compile_762*Compile_765;
      std::complex<double> Compile_1850 = \
      -0.0009765625*lam4H3dUS*Compile_340*Compile_418*Compile_494*Compile_747*Compile_762*Compile_765;
      std::complex<double> Compile_1851 = \
      0.0009765625*lam5H3dUS*Compile_340*Compile_418*Compile_494*Compile_747*Compile_762*Compile_765;
      std::complex<double> Compile_1823 = \
      0.0625*Compile_336*Compile_747*Compile_762*Compile_765;
      std::complex<double> Compile_1824 = 0. + Compile_1822 + Compile_1823;
      std::complex<double> Compile_1853 = \
      -1.*lam5H3dUS*Compile_527*Compile_807;
      std::complex<double> Compile_1862 = cos(Compile_806);
      std::complex<double> Compile_1318 = 0.5*Compile_816;
      std::complex<double> Compile_1325 = 0.5*Compile_949;
      std::complex<double> Compile_1863 = phi2*Compile_1862*Compile_768;
      std::complex<double> Compile_1864 = 0. + lam3H3dUS + lam4H3dUS + \
      Compile_1010 + Compile_410;
      std::complex<double> Compile_1865 = \
      -1.*phi1*Compile_1864*Compile_807;
      std::complex<double> Compile_1866 = Compile_1863 + Compile_1865;
      std::complex<double> Compile_1867 = \
      0.125*Compile_1866*Compile_417*Compile_439*Compile_555;
      std::complex<double> Compile_1868 = \
      -1.*phi1*Compile_1862*Compile_768;
      std::complex<double> Compile_1869 = 0. + lam3H3dUS + lam4H3dUS + \
      Compile_1011 + Compile_410;
      std::complex<double> Compile_1870 = \
      -1.*phi2*Compile_1869*Compile_807;
      std::complex<double> Compile_1871 = Compile_1868 + Compile_1870;
      std::complex<double> Compile_1872 = Compile_1871*Compile_555;
      std::complex<double> Compile_1873 = Compile_1867 + Compile_1872;
      std::complex<double> Compile_1874 = pow(Compile_1873,2);
      std::complex<double> Compile_1881 = -1.*Compile_1866*Compile_555;
      std::complex<double> Compile_1882 = \
      0.125*Compile_1871*Compile_417*Compile_439*Compile_555;
      std::complex<double> Compile_1883 = Compile_1881 + Compile_1882;
      std::complex<double> Compile_1884 = pow(Compile_1883,2);
      std::complex<double> Compile_1898 = Compile_1317 + Compile_1318 + \
      Compile_1422;
      std::complex<double> Compile_1913 = 0.25*Compile_815;
      std::complex<double> Compile_1893 = 1.*Compile_673;
      std::complex<double> Compile_1894 = -Compile_800;
      std::complex<double> Compile_1895 = Compile_1893 + Compile_1894;
      std::complex<double> Compile_1931 = Compile_1312 + Compile_1318 + \
      Compile_1422;
      std::complex<double> Compile_1917 = -1.*Compile_765*Compile_794;
      std::complex<double> Compile_1918 = \
      -0.0625*Compile_747*Compile_762*Compile_765*Compile_797;
      std::complex<double> Compile_1919 = Compile_1917 + Compile_1918 + \
      Compile_942;
      std::complex<double> Compile_1920 = 0.25*Compile_1919;
      std::complex<double> Compile_1929 = cos(Compile_1895);
      std::complex<double> Compile_1930 = pow(Compile_1929,2);
      std::complex<double> Compile_1963 = Compile_1317 + Compile_1325 + \
      Compile_1422;
      std::complex<double> Compile_1978 = 0.25*Compile_948;
      std::complex<double> Compile_1896 = sin(Compile_1895);
      std::complex<double> Compile_1897 = pow(Compile_1896,2);
      std::complex<double> Compile_1957 = -1.*Compile_765*Compile_797;
      std::complex<double> Compile_1958 = \
      -0.0625*Compile_747*Compile_762*Compile_765*Compile_794;
      std::complex<double> Compile_1959 = Compile_1957 + Compile_1958 + \
      Compile_808;
      std::complex<double> Compile_1960 = 0.25*Compile_1959;
      std::complex<double> Compile_1990 = Compile_1312 + Compile_1325 + \
      Compile_1422;
      std::complex<double> Compile_2023 = -0.5*Compile_949;
      std::complex<double> Compile_2022 = -0.5*Compile_269*Compile_273;
      std::complex<double> Compile_2017 = Compile_1325 + Compile_274 + \
      Compile_695;
      std::complex<double> Compile_2014 = Compile_525 + Compile_800;
      std::complex<double> Compile_2015 = cos(Compile_2014);
      std::complex<double> Compile_2016 = pow(Compile_2015,2);
      std::complex<double> Compile_2051 = -0.5*Compile_816;
      std::complex<double> Compile_2046 = Compile_1318 + Compile_274 + \
      Compile_735;
      std::complex<double> Compile_2074 = Compile_1325 + Compile_274 + \
      Compile_735;
      std::complex<double> Compile_2039 = \
      0.009947183943243459*Compile_61*Compile_64*Compile_949;
      std::complex<double> Compile_2102 = Compile_1318 + Compile_274 + \
      Compile_695;
      std::complex<double> Compile_2067 = \
      0.009947183943243459*Compile_61*Compile_64*Compile_816;
      std::complex<double> Compile_2099 = sin(Compile_2014);
      std::complex<double> Compile_2100 = pow(Compile_2099,2);
      return 0. + \
      0.026525823848649224*(1.4142135623730951*pow(Compile_27,1.5) + \
      1.4142135623730951*pow(Compile_31,1.5)) + 0.5*Compile_12*Compile_44 + \
      0.25*lam1H3dUS*Compile_49 + phi1*phi2*Compile_51 + \
      0.5*Compile_13*Compile_53 + 0.25*lam3H3dUS*Compile_44*Compile_53 + \
      0.25*lam4H3dUS*Compile_44*Compile_53 + \
      0.25*lam5H3dUS*Compile_44*Compile_53 - \
      0.003315727981081153*pow(Compile_534,1.5) - \
      0.003315727981081153*pow(Compile_546,1.5) + 0.25*lam2H3dUS*Compile_58 \
      + 0.006332573977646111*Compile_60*Compile_61 + \
      0.00019789293680144097*(1. + 2.*Compile_129)*Compile_60*Compile_61 - \
      0.019894367886486918*pow(Compile_65,1.5) - \
      0.006631455962162306*pow(Compile_692,1.5) - \
      0.006631455962162306*pow(Compile_723,1.5) - \
      0.0031662869888230555*(-0.5*(0. + Compile_1283 + \
      0.0625*Compile_1279*Compile_1284*Compile_1288*Compile_336 + \
      Compile_1279*Compile_450)*Compile_692 - 0.5*(0. + Compile_1283 + \
      Compile_1279*Compile_336 + \
      0.0625*Compile_1279*Compile_1284*Compile_1288*Compile_450)*Compile_723) - 0.009498860966469166*(-0.25*Compile_1270*Compile_64*Compile_693 - 0.25*Compile_1270*Compile_64*Compile_724) + \
      0.00039578587360288194*Compile_60*Compile_61*(1. + 2.*Compile_75) - \
      0.003315727981081153*pow(Compile_815,1.5) - \
      0.003315727981081153*pow(Compile_948,1.5) - \
      0.0007915717472057639*(-0.75*(0. + Compile_1710 + \
      Compile_1708*Compile_336 + \
      0.000244140625*Compile_1708*Compile_1711*Compile_1712*Compile_450)*\
      Compile_534 - 0.75*(0. + Compile_1710 + \
      0.000244140625*Compile_1708*Compile_1711*Compile_1712*Compile_336 + \
      Compile_1708*Compile_450)*Compile_546 - 0.75*(0. + Compile_1722 + \
      0.00390625*Compile_1720*Compile_1723*Compile_1724*Compile_336 + \
      Compile_1720*Compile_450)*Compile_815 - 0.75*(0. + Compile_1722 + \
      Compile_1720*Compile_336 + \
      0.00390625*Compile_1720*Compile_1723*Compile_1724*Compile_450)*\
      Compile_948) - \
      0.004749430483234583*(-0.25*Compile_1270*Compile_64*Compile_694 - \
      0.25*Compile_1270*Compile_64*Compile_734 - \
      0.25*Compile_1270*Compile_64*Compile_816 - \
      0.25*Compile_1270*Compile_64*Compile_949) - \
      0.0031662869888230555*(0.125*(0. + Compile_1764 + Compile_1766 + \
      Compile_1767 - 1.*Compile_1750*Compile_494 - \
      0.015625*Compile_340*Compile_418*Compile_494*Compile_557*Compile_608)*\
      Compile_693*Compile_694 + 0.25*(-0.5*(0. + Compile_1743 + \
      Compile_1744)*Compile_494 - \
      0.0078125*Compile_1741*Compile_340*Compile_418*Compile_494 + \
      0.125*Compile_417*Compile_439*Compile_494*(0. + \
      0.5*Compile_600*Compile_608*Compile_611*Compile_612))*Compile_694*\
      Compile_724 + 0.125*(-1.*(0. + Compile_1744 + \
      Compile_1752)*Compile_494 - \
      0.015625*Compile_1750*Compile_340*Compile_418*Compile_494 - \
      0.25*Compile_417*Compile_439*Compile_494*(0. - 1.*Compile_1755*(0. + \
      Compile_1755*Compile_559)*Compile_600*Compile_612))*Compile_693*\
      Compile_734 + 0.125*(0. + Compile_1764 + Compile_1766 + Compile_1767 \
      - 1.*Compile_1741*Compile_494 - \
      0.015625*Compile_340*Compile_418*Compile_494*Compile_608*Compile_621)*\
      Compile_724*Compile_734 + 0.125*Compile_724*(0. + Compile_1802 + \
      Compile_1804 + Compile_1805 - 1.*Compile_1789*Compile_765 - \
      0.0625*Compile_608*Compile_621*Compile_747*Compile_762*Compile_765)*\
      Compile_816 + 0.125*Compile_693*(-1.*(0. + Compile_1752 + \
      0.25*Compile_570*Compile_601*Compile_608*Compile_614)*Compile_765 - \
      0.0625*Compile_1775*Compile_747*Compile_762*Compile_765 + 0.5*(0. - \
      1.*lam4H3dUS*Compile_683 - \
      1.*lam5H3dUS*Compile_683)*Compile_761*Compile_765*Compile_769)*\
      Compile_816 + 0.125*Compile_693*(0. + Compile_1802 + Compile_1804 + \
      Compile_1805 - 1.*Compile_1775*Compile_765 - \
      0.0625*Compile_557*Compile_608*Compile_747*Compile_762*Compile_765)*\
      Compile_949 + 0.125*Compile_724*(-1.*(Compile_1743 - \
      0.5*Compile_1755*Compile_600*Compile_612*(0. - \
      0.5*Compile_1755*Compile_600*Compile_612*Compile_614))*Compile_765 - \
      0.0625*Compile_1789*Compile_747*Compile_762*Compile_765 + 0.25*(0. - \
      2.*Compile_559*Compile_683)*Compile_761*Compile_765*Compile_769)*\
      Compile_949) - 1.*(Compile_1245 + Compile_1246 + \
      0.0015831434944115277*Compile_1032*Compile_64 - \
      0.0015831434944115277*Compile_1041*Compile_64 + \
      0.0015831434944115277*Compile_1053*(0.5*(2.*Compile_1008*Compile_1009 \
      + Compile_1061 + Compile_1062) + \
      2.*Compile_1063*Compile_1073)*Compile_64 - \
      0.0015831434944115277*Compile_1053*Compile_1075*Compile_64 + \
      0.5*(Compile_1245 + Compile_1246 - \
      0.0031662869888230555*Compile_1053*Compile_1075*Compile_64) - \
      0.0031662869888230555*(-0.375*Compile_27*Compile_64 - \
      0.5*Compile_1029*Compile_27*Compile_64 - 0.375*Compile_31*Compile_64 \
      - 0.5*Compile_1039*Compile_31*Compile_64) - \
      0.0007915717472057639*Compile_1008*Compile_1009*(Compile_1010 + \
      Compile_1011 + Compile_1012 - 2.*Compile_1015*Compile_1020 + \
      Compile_585 + Compile_767) - \
      0.0031662869888230555*(-1.*Compile_31*Compile_979 - \
      1.*Compile_27*Compile_983) - \
      0.0031662869888230555*(0.25*Compile_31*(-1.*Compile_1088*Compile_1165 \
      - 1.*Compile_1016*Compile_1088*Compile_1140*Compile_559*Compile_969 - \
      0.25*Compile_1088*Compile_1159*Compile_971*Compile_972) + \
      0.25*Compile_27*(-1.*Compile_1088*Compile_1174 + \
      Compile_1016*Compile_1088*(-1.*lam4H3dUS*Compile_1140 - \
      1.*lam5H3dUS*Compile_1140)*Compile_969 - \
      0.25*Compile_1088*Compile_1171*Compile_971*Compile_972) + \
      0.25*Compile_1008*Compile_1009*(-1.*Compile_1088*Compile_1171 + \
      Compile_1189 + Compile_1194 + Compile_1195 - \
      1.*lam1H3dUS*Compile_971*Compile_972*Compile_975) + \
      0.25*Compile_1008*Compile_1009*(-1.*Compile_1088*Compile_1159 + \
      Compile_1189 + Compile_1194 + Compile_1195 - \
      1.*lam2H3dUS*Compile_971*Compile_972*Compile_975) + \
      0.25*Compile_27*(Compile_1147 - 1.*Compile_1087*Compile_1171 + \
      Compile_1215 + Compile_1217 - \
      1.*lam1H3dUS*Compile_1087*Compile_1088*Compile_971*Compile_992) + \
      0.25*Compile_31*(Compile_1147 - 1.*Compile_1087*Compile_1159 + \
      Compile_1215 + Compile_1217 - \
      1.*lam2H3dUS*Compile_1087*Compile_1088*Compile_971*Compile_992) + \
      0.5*Compile_1008*Compile_1009*(-0.5*Compile_1087*Compile_1165 + \
      0.5*Compile_1087*Compile_1088*Compile_559*Compile_969*Compile_971*\
      Compile_990 - \
      0.125*Compile_1087*Compile_1159*Compile_971*Compile_992) + \
      0.25*Compile_1008*Compile_1009*(-1.*Compile_1087*Compile_1174 + \
      Compile_1087*Compile_1088*Compile_559*Compile_969*Compile_971*Compile_990 - 0.25*Compile_1087*Compile_1171*Compile_971*Compile_992)) - \
      0.0007915717472057639*(-1.5*Compile_31*Compile_979 - \
      1.5*Compile_27*Compile_983 - 1.5*Compile_27*(2.*lam2H3dUS*Compile_995 \
      + 0.125*lam1H3dUS*Compile_968*Compile_991*Compile_995 + Compile_998) \
      - 1.5*Compile_31*(2.*lam1H3dUS*Compile_995 + \
      0.125*lam2H3dUS*Compile_968*Compile_991*Compile_995 + Compile_998)) - \
      0.0015831434944115277*(0.0625*Compile_1008*Compile_1009*(Compile_1012 \
      - 6.*Compile_1015*Compile_1020 + Compile_1119 + Compile_1121 + \
      Compile_399 + Compile_585) + \
      0.5*Compile_27*(-1.*Compile_1087*Compile_1091 + Compile_1135 + \
      Compile_1136 + Compile_1138 + Compile_1147 - \
      0.5*lam1H3dUS*Compile_1087*Compile_1088*Compile_971*Compile_992) + \
      0.5*Compile_31*(-1.*Compile_1087*Compile_1112 + Compile_1135 + \
      Compile_1136 + Compile_1138 + Compile_1147 - \
      0.5*lam2H3dUS*Compile_1087*Compile_1088*Compile_971*Compile_992) + \
      0.5*Compile_1008*Compile_1009*(-1.*Compile_1087*(2.*lam1H3dUS*Compile_1088 + Compile_1101) + Compile_1104 - \
      0.25*Compile_1087*Compile_1091*Compile_971*Compile_992) + \
      0.5*Compile_1008*Compile_1009*(-1.*Compile_1087*(2.*lam2H3dUS*Compile_1088 + Compile_1101) + Compile_1104 - \
      0.25*Compile_1087*Compile_1112*Compile_971*Compile_992) + \
      0.0625*Compile_1008*Compile_1009*(Compile_1012 + Compile_1119 + \
      Compile_1121 + Compile_399 + Compile_585 - \
      6.*Compile_1015*cos(4.*Compile_1050)))) - \
      0.00039578587360288194*Compile_693*Compile_724*(0. + Compile_1010 + \
      Compile_1011 + Compile_1012 + Compile_585 + Compile_767 - \
      1.*Compile_1261*cos(4.*Compile_673)) - \
      0.0015831434944115277*(0.25*Compile_694*(0. + Compile_1849 + \
      Compile_1850 + Compile_1851 + Compile_1853 - \
      1.*Compile_1824*Compile_494 - \
      0.015625*Compile_340*Compile_418*Compile_450*Compile_494*Compile_765)*\
      Compile_816 + \
      0.25*Compile_734*(-0.015625*Compile_1824*Compile_340*Compile_418*\
      Compile_494 - 1.*Compile_494*(0. + Compile_1827 + \
      Compile_450*Compile_765) - \
      0.25*Compile_417*Compile_439*Compile_494*(0. - \
      0.25*Compile_761*Compile_765*Compile_768*Compile_769))*Compile_816 + \
      0.25*Compile_734*(0. + Compile_1849 + Compile_1850 + Compile_1851 + \
      Compile_1853 - 1.*Compile_1836*Compile_494 - \
      0.015625*Compile_336*Compile_340*Compile_418*Compile_494*Compile_765)*\
      Compile_949 + \
      0.25*Compile_694*(-0.015625*Compile_1836*Compile_340*Compile_418*\
      Compile_494 - 1.*Compile_494*(0. + Compile_1827 + \
      Compile_336*Compile_765) + \
      0.25*Compile_417*Compile_439*Compile_494*(0. + \
      0.25*Compile_761*Compile_765*Compile_768*Compile_769))*Compile_949 + \
      0.03125*Compile_694*Compile_734*(0. + Compile_1012 + Compile_1119 + \
      Compile_1121 + Compile_399 + Compile_585 - \
      3.*Compile_1261*cos(4.*Compile_525)) + \
      0.03125*Compile_816*Compile_949*(0. + Compile_1012 + Compile_1119 + \
      Compile_1121 + Compile_399 + Compile_585 - \
      3.*Compile_1261*cos(4.*Compile_800))) - \
      1.*Compile_1368*Compile_1897*(0.0015831434944115277*Compile_1270*(\
      Compile_1450 + Compile_1451 + Compile_1920)*Compile_724 - \
      0.039788735772973836*(Compile_1461 - \
      0.039788735772973836*Compile_1270*(Compile_1445 + Compile_1450 + \
      Compile_1913))*Compile_816 - 0.0031662869888230555*pow(Compile_1445 + \
      Compile_1913,2)*(1. + 2.*log(mu3US/(Compile_1317 + Compile_1318))) + \
      0.006332573977646111*(Compile_1317 + Compile_1318 + \
      Compile_1429)*(Compile_1318 + Compile_1422 + \
      Compile_1431)*(Compile_1318 + Compile_1429 + \
      Compile_1431)*Compile_1898*(0.5 + log(mu3US/Compile_1898))) - \
      1.*Compile_1368*Compile_1930*(0.0015831434944115277*Compile_1270*(\
      Compile_1450 + Compile_1474 + Compile_1920)*Compile_693 - \
      0.039788735772973836*(Compile_1494 - \
      0.039788735772973836*Compile_1270*(Compile_1450 + Compile_1498 + \
      Compile_1913))*Compile_816 - 0.0031662869888230555*pow(Compile_1498 + \
      Compile_1913,2)*(1. + 2.*log(mu3US/(Compile_1312 + Compile_1318))) + \
      0.006332573977646111*(Compile_1312 + Compile_1318 + \
      Compile_1429)*(Compile_1318 + Compile_1422 + \
      Compile_1489)*(Compile_1318 + Compile_1429 + \
      Compile_1489)*Compile_1931*(0.5 + log(mu3US/Compile_1931))) - \
      1.*Compile_1368*Compile_1930*(0.0015831434944115277*Compile_1270*(\
      Compile_1450 + Compile_1451 + Compile_1960)*Compile_724 - \
      0.039788735772973836*(Compile_1461 - \
      0.039788735772973836*Compile_1270*(Compile_1445 + Compile_1450 + \
      Compile_1978))*Compile_949 - 0.0031662869888230555*pow(Compile_1445 + \
      Compile_1978,2)*(1. + 2.*log(mu3US/(Compile_1317 + Compile_1325))) + \
      0.006332573977646111*(Compile_1317 + Compile_1325 + \
      Compile_1429)*(Compile_1325 + Compile_1422 + \
      Compile_1431)*(Compile_1325 + Compile_1429 + \
      Compile_1431)*Compile_1963*(0.5 + log(mu3US/Compile_1963))) - \
      1.*Compile_1368*Compile_1897*(0.0015831434944115277*Compile_1270*(\
      Compile_1450 + Compile_1474 + Compile_1960)*Compile_693 - \
      0.039788735772973836*(Compile_1494 - \
      0.039788735772973836*Compile_1270*(Compile_1450 + Compile_1498 + \
      Compile_1978))*Compile_949 - 0.0031662869888230555*pow(Compile_1498 + \
      Compile_1978,2)*(1. + 2.*log(mu3US/(Compile_1312 + Compile_1325))) + \
      0.006332573977646111*(Compile_1312 + Compile_1325 + \
      Compile_1429)*(Compile_1325 + Compile_1422 + \
      Compile_1489)*(Compile_1325 + Compile_1429 + \
      Compile_1489)*Compile_1990*(0.5 + log(mu3US/Compile_1990))) - \
      0.5*Compile_64*(0.009103075092866285*Compile_61*Compile_64 + \
      0.00039578587360288194*(0.5 + Compile_129)*Compile_61*Compile_64 + \
      0.0007915717472057639*Compile_61*Compile_64*(0.5 + Compile_75) + \
      0.006332573977646111*((0.125*(9.*Compile_44*Compile_60 + \
      9.*Compile_53*Compile_60))/pow(g23dUS,2) + \
      (0.0625*(-63.*Compile_260*Compile_44 - \
      63.*Compile_260*Compile_53)*(0.5 + log(mu3US/(Compile_270 + \
      Compile_274))))/pow(g23dUS,4))) - \
      0.00026385724906858796*(9.*Compile_456 + 9.*Compile_482 + \
      18.*Compile_456*log(0.6666666666666666*mu3US*Compile_536) + \
      18.*Compile_482*log(0.6666666666666666*mu3US*Compile_547)) - \
      1.*Compile_1368*Compile_1670*(0.0015831434944115277*Compile_1270*(\
      Compile_1450 + Compile_1457 + Compile_1474)*Compile_693 - \
      0.039788735772973836*(Compile_1494 - \
      0.039788735772973836*Compile_1270*(Compile_1446 + Compile_1450 + \
      Compile_1498))*Compile_694 + \
      0.006332573977646111*Compile_1644*(Compile_1312 + Compile_1429 + \
      Compile_695)*(Compile_1422 + Compile_1489 + \
      Compile_695)*(Compile_1429 + Compile_1489 + Compile_695)*(0.5 + \
      log(mu3US/Compile_1644)) - 0.0031662869888230555*pow(Compile_1457 + \
      Compile_1474,2)*(1. + 2.*log(mu3US/(Compile_1312 + Compile_695)))) - \
      1.*Compile_1368*Compile_1421*(-0.039788735772973836*(-0.039788735772973836*Compile_1270*(Compile_1445 + Compile_1446 + \
      Compile_1450) + Compile_1461)*Compile_694 + \
      0.0015831434944115277*Compile_1270*(Compile_1450 + Compile_1451 + \
      Compile_1457)*Compile_724 + \
      0.006332573977646111*Compile_1424*(Compile_1317 + Compile_1429 + \
      Compile_695)*(Compile_1422 + Compile_1431 + \
      Compile_695)*(Compile_1429 + Compile_1431 + Compile_695)*(0.5 + \
      log(mu3US/Compile_1424)) - 0.0031662869888230555*pow(Compile_1445 + \
      Compile_1446,2)*(1. + 2.*log(mu3US/(Compile_1317 + Compile_695)))) + \
      0.5*(-1.*Compile_1368*(0.00039578587360288194*Compile_1270*Compile_61*\
      Compile_64*Compile_693 + \
      0.00039578587360288194*Compile_61*Compile_64*(Compile_1270 - \
      1.*Compile_693)*Compile_693 + \
      0.00019789293680144097*Compile_61*Compile_64*(Compile_1385 - \
      2.*Compile_693)*(Compile_1385 + 2.*Compile_693)*(1. + \
      2.*log(mu3US/(Compile_274 + Compile_693)))) - \
      1.*Compile_1368*(0.00039578587360288194*Compile_1270*Compile_61*\
      Compile_64*Compile_724 + \
      0.00039578587360288194*Compile_61*Compile_64*(Compile_1270 - \
      1.*Compile_724)*Compile_724 + \
      0.00019789293680144097*Compile_61*Compile_64*(Compile_1385 - \
      2.*Compile_724)*(Compile_1385 + 2.*Compile_724)*(1. + \
      2.*log(mu3US/(Compile_274 + Compile_724))))) - \
      1.*Compile_1368*Compile_1421*(0.0015831434944115277*Compile_1270*(\
      Compile_1450 + Compile_1474 + Compile_1478)*Compile_693 - \
      0.039788735772973836*(Compile_1494 - \
      0.039788735772973836*Compile_1270*(Compile_1450 + Compile_1498 + \
      Compile_1499))*Compile_734 + \
      0.006332573977646111*Compile_1484*(Compile_1312 + Compile_1429 + \
      Compile_735)*(Compile_1422 + Compile_1489 + \
      Compile_735)*(Compile_1429 + Compile_1489 + Compile_735)*(0.5 + \
      log(mu3US/Compile_1484)) - 0.0031662869888230555*pow(Compile_1474 + \
      Compile_1478,2)*(1. + 2.*log(mu3US/(Compile_1312 + Compile_735)))) - \
      1.*Compile_1368*Compile_1670*(0.0015831434944115277*Compile_1270*(\
      Compile_1450 + Compile_1451 + Compile_1478)*Compile_724 - \
      0.039788735772973836*(Compile_1461 - \
      0.039788735772973836*Compile_1270*(Compile_1445 + Compile_1450 + \
      Compile_1499))*Compile_734 + \
      0.006332573977646111*Compile_1674*(Compile_1317 + Compile_1429 + \
      Compile_735)*(Compile_1422 + Compile_1431 + \
      Compile_735)*(Compile_1429 + Compile_1431 + Compile_735)*(0.5 + \
      log(mu3US/Compile_1674)) - 0.0031662869888230555*pow(Compile_1445 + \
      Compile_1499,2)*(1. + 2.*log(mu3US/(Compile_1317 + Compile_735)))) - \
      0.0031662869888230555*(0.0625*Compile_1301*Compile_1305 + \
      0.0625*Compile_1301*Compile_1310 + 0.0625*Compile_1349 + \
      0.0625*Compile_1359 + \
      0.125*Compile_1301*Compile_1305*log(mu3US/(Compile_1312 + \
      Compile_1317 + Compile_1318)) + \
      0.125*Compile_1301*Compile_1310*log(mu3US/(Compile_1312 + \
      Compile_1317 + Compile_1325)) + \
      0.125*Compile_1349*log(mu3US/(Compile_1312 + Compile_1317 + \
      Compile_695)) + 0.125*Compile_1359*log(mu3US/(Compile_1312 + \
      Compile_1317 + Compile_735))) + \
      0.5*(-1.*Compile_1368*Compile_2100*(-0.039788735772973836*(-0.039788735772973836*Compile_1270*(Compile_1446 + Compile_1450 + \
      Compile_1920) + Compile_2067)*Compile_694 + \
      0.0015831434944115277*Compile_1270*(Compile_1450 + Compile_1457 + \
      Compile_1913)*Compile_816 + \
      0.006332573977646111*Compile_2102*(Compile_1318 + Compile_2022 + \
      Compile_695)*(Compile_2022 + Compile_2051 + \
      Compile_695)*(Compile_2051 + Compile_274 + Compile_695)*(0.5 + \
      log(mu3US/Compile_2102)) - 0.0031662869888230555*pow(Compile_1457 + \
      Compile_1913,2)*(1. + 2.*log(mu3US/(Compile_1318 + Compile_695)))) - \
      1.*Compile_1368*Compile_2016*(-0.039788735772973836*(-0.039788735772973836*Compile_1270*(Compile_1446 + Compile_1450 + \
      Compile_1960) + Compile_2039)*Compile_694 + \
      0.0015831434944115277*Compile_1270*(Compile_1450 + Compile_1457 + \
      Compile_1978)*Compile_949 + \
      0.006332573977646111*Compile_2017*(Compile_1325 + Compile_2022 + \
      Compile_695)*(Compile_2022 + Compile_2023 + \
      Compile_695)*(Compile_2023 + Compile_274 + Compile_695)*(0.5 + \
      log(mu3US/Compile_2017)) - 0.0031662869888230555*pow(Compile_1457 + \
      Compile_1978,2)*(1. + 2.*log(mu3US/(Compile_1325 + Compile_695)))) - \
      1.*Compile_1368*Compile_2016*(-0.039788735772973836*(-0.039788735772973836*Compile_1270*(Compile_1450 + Compile_1499 + \
      Compile_1920) + Compile_2067)*Compile_734 + \
      0.0015831434944115277*Compile_1270*(Compile_1450 + Compile_1478 + \
      Compile_1913)*Compile_816 + \
      0.006332573977646111*Compile_2046*(Compile_1318 + Compile_2022 + \
      Compile_735)*(Compile_2022 + Compile_2051 + \
      Compile_735)*(Compile_2051 + Compile_274 + Compile_735)*(0.5 + \
      log(mu3US/Compile_2046)) - 0.0031662869888230555*pow(Compile_1478 + \
      Compile_1913,2)*(1. + 2.*log(mu3US/(Compile_1318 + Compile_735)))) - \
      1.*Compile_1368*Compile_2100*(-0.039788735772973836*(-0.039788735772973836*Compile_1270*(Compile_1450 + Compile_1499 + \
      Compile_1960) + Compile_2039)*Compile_734 + \
      0.0015831434944115277*Compile_1270*(Compile_1450 + Compile_1478 + \
      Compile_1978)*Compile_949 + \
      0.006332573977646111*Compile_2074*(Compile_1325 + Compile_2022 + \
      Compile_735)*(Compile_2022 + Compile_2023 + \
      Compile_735)*(Compile_2023 + Compile_274 + Compile_735)*(0.5 + \
      log(mu3US/Compile_2074)) - 0.0031662869888230555*pow(Compile_1478 + \
      Compile_1978,2)*(1. + 2.*log(mu3US/(Compile_1325 + Compile_735))))) - \
      0.0015831434944115277*(0.25*Compile_1874 + 0.25*Compile_1884 + \
      0.5*Compile_1874*log(mu3US/(Compile_1318 + Compile_1325 + \
      Compile_695)) + 0.5*Compile_1884*log(mu3US/(Compile_1318 + \
      Compile_1325 + Compile_735))) + \
      0.5*(-1.*Compile_1506*Compile_1558*(Compile_1559 + Compile_1565 + \
      Compile_1591 + 0.012665147955292222*Compile_1573*(-0.5 - \
      1.*log(mu3US/(Compile_1422 + Compile_695))) + \
      0.006332573977646111*(Compile_1537 + (Compile_1543 + Compile_1544 + \
      Compile_1545 + 0.0625*Compile_1564 - \
      0.25*Compile_44*Compile_534*Compile_64 - \
      0.25*Compile_53*Compile_534*Compile_64)*(0.5 + log(mu3US/(Compile_270 \
      + Compile_695))))) - 1.*Compile_1506*Compile_1510*(Compile_1511 + \
      Compile_1517 + Compile_1532 + 0.012665147955292222*Compile_1525*(-0.5 \
      - 1.*log(mu3US/(Compile_1422 + Compile_735))) + \
      0.006332573977646111*(Compile_1537 + (0.0625*Compile_1516 + \
      Compile_1543 + Compile_1544 + Compile_1545 - \
      0.25*Compile_44*Compile_546*Compile_64 - \
      0.25*Compile_53*Compile_546*Compile_64)*(0.5 + log(mu3US/(Compile_270 \
      + Compile_735)))))) + \
      0.25*(-1.*Compile_1506*Compile_1558*(Compile_1559 + Compile_1565 + \
      Compile_1591 + 0.006332573977646111*(Compile_1604 + (Compile_1573 + \
      Compile_1610 + 0.25*(0.5*Compile_1454 + \
      Compile_1611)*Compile_61*Compile_64)*(0.5 + log(mu3US/(Compile_1385 + \
      Compile_695)))) + 0.012665147955292222*Compile_1573*(-0.5 - \
      1.*log(mu3US/(Compile_274 + Compile_695)))) - \
      1.*Compile_1506*Compile_1510*(Compile_1511 + Compile_1517 + \
      Compile_1532 + 0.006332573977646111*(Compile_1604 + (Compile_1525 + \
      Compile_1610 + 0.25*(0.5*Compile_1477 + \
      Compile_1611)*Compile_61*Compile_64)*(0.5 + log(mu3US/(Compile_1385 + \
      Compile_735)))) + 0.012665147955292222*Compile_1525*(-0.5 - \
      1.*log(mu3US/(Compile_274 + Compile_735))))) - \
      0.0015831434944115277*(0.25*Compile_628 + 0.25*Compile_633 + \
      0.25*Compile_651 + 0.25*Compile_657 + 0.5*Compile_628*Compile_700 + \
      0.5*Compile_633*Compile_700 + 0.5*Compile_651*log(mu3US/(Compile_695 \
      + Compile_724)) + 0.5*Compile_657*log(mu3US/(Compile_724 + \
      Compile_735))) + \
      0.25*(-0.0031662869888230555*pow(0.125*Compile_417*Compile_439*\
      Compile_555*Compile_773 - 1.*Compile_555*Compile_786,2)*(1. + \
      2.*log(mu3US/(Compile_695 + Compile_816))) - \
      0.0031662869888230555*(Compile_843 + Compile_875 + Compile_903 + \
      Compile_907 + Compile_922 + 2.*Compile_843*log(mu3US/(Compile_695 + \
      Compile_734)) + 2.*Compile_875*log(mu3US/(Compile_694 + Compile_735)) \
      + 2.*Compile_922*log(mu3US/(Compile_735 + Compile_816)) + \
      2.*Compile_903*log(mu3US/(Compile_695 + Compile_949)) + \
      2.*Compile_907*log(mu3US/(Compile_735 + Compile_949))));
    }


  };