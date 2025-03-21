/**
  Z2 real scalar singlet extension of
  the Standard Model 
  in MSbar scheme.
  See arXiv:2208.01319  [hep-ph] for details
  
  
*/

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <random>

#include "models/xSM_MSbar.hpp"
#include "phasetracer.hpp" 

std::string toString(std::vector<double> in, std::vector<double> out, std::vector<double> flags) {
  std::stringstream data_str;
  for (auto i : in    ) data_str << i << "\t";
  for (auto i : out   ) data_str << i << "\t";
  for (auto i : flags ) data_str << i << "\t";
  return data_str.str();;
}

int main(int argc, char* argv[]) {

  std::ofstream output_file;  
  output_file.open("output.txt");

  bool debug_mode = false;
  double ms, lambda_s, lambda_hs;
  double Q, xi, daisy_flag;
  bool use_1L_EWSB_in_0L_mass;
  bool use_Goldstone_resum = true;
  bool tree_level_tadpoles = false;
  bool use_covariant_gauge = false;
  std::vector<double> SM_parameters ={};
  if ( argc == 1 ) {
    debug_mode = true;
    // Compare with run_ScalarSingletZ2DMMhInput_withSingletVEVinPT
    ms = 84.6733668341708;
    lambda_s =  0.1;
    lambda_hs = 0.3;
    Q = 86.5;
    xi = 10.;
    daisy_flag = 2;
    use_1L_EWSB_in_0L_mass = false;
    use_Goldstone_resum = true;
    use_covariant_gauge = false;
    
  } else if ( argc >= 9 ) {
    ms = atof(argv[1]);
    lambda_s = atof(argv[2]);
    lambda_hs = atof(argv[3]);
    Q = atof(argv[4]);
    xi = atof(argv[5]);

    daisy_flag = atoi(argv[6]);
    use_1L_EWSB_in_0L_mass = atoi(argv[7]);
    use_Goldstone_resum = atoi(argv[8]);
    if ( argc > 9 ){
      // default
      // tree_level_tadpoles = false
      // use_covariant_gauge = true
      if ( atoi(argv[9]) == 1 ){
        tree_level_tadpoles = true;
      } else if ( atoi(argv[9]) == 2 ){
        use_covariant_gauge = true;
      } else if ( atoi(argv[9]) == 3 ){
        use_covariant_gauge = true;
        tree_level_tadpoles = true;
      } 
    }
  } else {
    std::cout << "Use ./run_xSM_MSbar ms lambda_s lambda_hs Q xi daisy_flag use_1L_EWSB_in_0L_mass use_Goldstone_resum" << std::endl;
    return 0;
  }

  std::vector<double> in ={ms, lambda_s, lambda_hs};
  std::vector<double> flags ={Q, xi, daisy_flag, (float)use_1L_EWSB_in_0L_mass, (float)use_Goldstone_resum};

  if (debug_mode){
    LOGGER(debug);
    std::cout << "ms = " << ms << std::endl
              << "lambda_s = " << lambda_s << std::endl
              << "lambda_hs = " << lambda_hs << std::endl
              << "Q = " << Q << std::endl
              << "xi = " << xi << std::endl
              << "daisy_term = " << ( daisy_flag == 0  ? "None" : ( daisy_flag == 1 ? "Parwani" : "ArnoldEspinosa")) << std::endl
              << "use 1-level ewsb in tree-level masses = " << use_1L_EWSB_in_0L_mass << std::endl
              << "use Goldstone resum = " << use_Goldstone_resum << std::endl;

  } else {
    LOGGER(fatal);
  }
  
  // Construct our model
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, use_covariant_gauge, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, tree_level_tadpoles, SM_parameters);
  if (debug_mode) std::cout << "1-L EWSB iteration converged = " << model.iteration_converged << std::endl;

  if (not model.iteration_converged){
    std::cout << "ms = " << ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "found 0 transition!" << std::endl;
    std::vector<double> out = {-100, 0, 0, 0, 0, 0};
    output_file << toString(in, out, flags) << std::endl;
    return 0;
  }
  

  // Choose Daisy method 
  if (daisy_flag == 0){
    model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  } else if (daisy_flag == 1){
    model.set_daisy_method(EffectivePotential::DaisyMethod::Parwani);
  } else if (daisy_flag == 2){
    model.set_daisy_method(EffectivePotential::DaisyMethod::ArnoldEspinosa);
  } else {
      std::cout << "Wrong daisy flag" << std::endl;
  }


  if (debug_mode) {
    std::cout << std::setprecision(16);
    std::cout << std::endl;
    std::cout << "@ after applying the 1L EWSB condition" << std::endl;
    std::cout << "muH2       = "<< model.get_muh_sq() << std::endl;
    std::cout << "muS2       = "<< model.get_mus_sq() << std::endl;   
    std::cout << "lambda_h   = "<< model.get_lambda_h() << std::endl;      
    std::cout << "lambda_s   = "<< model.get_lambda_s() << std::endl;  
    std::cout << "lambda_hs  = "<< model.get_lambda_hs() << std::endl; 
    std::cout << "tree min   = "<< std::sqrt(-model.get_muh_sq() / model.get_lambda_h()) << std::endl; 
  
    
    std::cout << std::endl;
    std::cout << "@ EWSB VEV" << std::endl;
    Eigen::VectorXd test(2);
    test <<  0, 0;
    double Ttest = 0;
    
    auto mh_check =  model.get_scalar_masses_sq(test,0);
    auto mV_check = model.get_vector_masses_sq(test);
    auto mf_check = model.get_fermion_masses_sq(test);
//    std::cout << "mh1 = "<< std::sqrt(std::abs(mh_check[0])) << std::endl;
//    std::cout << "mh2 = "<< std::sqrt(std::abs(mh_check[1])) << std::endl;
//    std::cout << "mh3 = "<< std::sqrt(std::abs(mh_check[2])) << std::endl;
//    std::cout << "mh4 = "<< std::sqrt(std::abs(mh_check[3])) << std::endl;
//    std::cout << "mh5 = "<< std::sqrt(std::abs(mh_check[4])) << std::endl;
//    std::cout << "mh6 = "<< std::sqrt(std::abs(mh_check[5])) << std::endl;
//      return 0;
    std::cout << "MW = "<< std::sqrt(mV_check[0]) << std::endl;
    std::cout << "MZ = "<< std::sqrt(mV_check[1]) << std::endl;
    std::cout << "Mphoton = "<< std::sqrt(mV_check[2]) << std::endl;
    std::cout << "mt = " << std::sqrt(mf_check[0]) << std::endl;
    std::cout << "mb = " << std::sqrt(mf_check[1]) << std::endl;
    std::cout << "mtau = " << std::sqrt(mf_check[2]) << std::endl;
//    
    double Vtree = model.V0(test);
    double VCW = model.V1(test);
    double V1T = model.V1T(test, Ttest);
    double Vtot = model.V(test, Ttest);
    std::cout << "Vtree      = "<< Vtree << std::endl;
    std::cout << "VCW        = "<< VCW << std::endl;   
//    std::cout << "V1T(T=100) = "<< V1T << std::endl;      
//    std::cout << "V(T=100)   = "<< Vtot << std::endl;  
//    std::cout << std::endl;
    
//    std::cout << "Numerically derivatives of the full potential at EW VEV:" << std::endl;
//    auto d2Vdh2 = model.d2V_dx2(test,0);
//    std::cout << std::setprecision(16);
//    std::cout << "Sqrt[d^2V/dh^2] = "<< std::sqrt(abs(d2Vdh2(0,0))) << std::endl;
//    std::cout << "Sqrt[d^2V/ds^2] = "<< std::sqrt(abs(d2Vdh2(1,1))) << std::endl;



      std::string prefix;
      if (daisy_flag == 0){
        prefix = "nodaisy";
      } else if (daisy_flag == 1){
        prefix = "Parwani";
      } else if (daisy_flag == 2){
        prefix = "ArnoldEspinosa";
      }

      // scale
      Eigen::VectorXd x1(2), x2(2);
      x1 << 245, 0;
      x2 << 275, 0;
      PhaseTracer::potential_line_plotter(model, 0, x1 , x2, "no_RGE_0_1_"+std::to_string(Q));

//      Eigen::VectorXd x1(2), x2(2);
//      x1 << 0, 0;
//      x2 << 250, 0;
//      PhaseTracer::potential_line_plotter(model, 150, x1 , x2, "150_1_"+prefix);
//      x1 << 0, 0;
//      x2 << 0, 250;
//      PhaseTracer::potential_line_plotter(model, 150, x1 , x2, "150_2_"+prefix); 
//      
//      x1 << 0, 0;
//      x2 << 250, 0;
//      PhaseTracer::potential_line_plotter(model, 110.4339097933532, x1 , x2, "110_1_"+prefix);
//      x1 << 0, 0;
//      x2 << 0, 250;
//      PhaseTracer::potential_line_plotter(model, 110.4339097933532, x1 , x2, "110_2_"+prefix);
//      
//      
//      x1 << 0, 0;
//      x2 << 2, 0;
//      PhaseTracer::potential_line_plotter(model, 10, x1 , x2, "10_1_"+prefix);
//      x1 << 245, 0;
//      x2 << 247, 0;
//      PhaseTracer::potential_line_plotter(model, 10, x1 , x2, "10_2_"+prefix);

//      x1 << 0, 0;
//      x2 << 0, 2;
//      PhaseTracer::potential_line_plotter(model, 10, x1 , x2, "10_3_"+prefix);
//    
//      x1 << 0, 217;
//      x2 << 0, 219;
//      PhaseTracer::potential_line_plotter(model, 10, x1 , x2, "10_4_"+prefix);
    
//    PhaseTracer::potential_plotter(model, 254, "potential", -5., 5, 0.01, -5., 40., 0.1);
//    PhaseTracer::potential_plotter(model, 142.35, "potential", 0., 160, 0.2, -2., 160., 0.2);
//    return 0;
  }
      
  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);    
  pf.set_check_vacuum_at_high(false);
  pf.set_seed(1);
//  pf.set_check_hessian_singular(false);
//  pf.set_hessian_singular_rel_tol(1.e-6);
    
  try {
    pf.find_phases();
  } catch (...) {
    std::cout << "ms = " << ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "encounters bug!" << std::endl;
    std::vector<double> out = {-1, 0, 0, 0, 0, 0};
    output_file << toString(in, out, flags) << std::endl;
    return 0;
  }
    
  if (debug_mode) std::cout << pf;
  
  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  if (debug_mode) std::cout << tf;  
  
  auto t = tf.get_transitions();
  if (t.size()==0){
    std::cout << "ms = " << ms << ",\t"
              << "lambda_s = " << lambda_s << ",\t"
              << "lambda_hs = " << lambda_hs << "\t"
              << "found 0 transition!" << std::endl;
    std::vector<double> out = {-2, 0, 0, 0, 0, 0};
    output_file << toString(in, out, flags) << std::endl;
    return 0;
  }
    
  // Find the transition with largest gamma from (0,vs) -> (vh,0) 
  int jj = -1;
  double gamme_max = 0.;
  for (int i=0; i<t.size(); i++) {
    double gamma = t[i].gamma;
    if (gamme_max < gamma and abs(t[i].true_vacuum[1])<1. and abs(t[i].false_vacuum[0])<1.){
      jj = i;
      gamme_max = gamma;
    }
  }
  
  if (jj<0) {
    std::vector<double> out = {-3, 0, 0, 0, 0, 0};
    output_file << toString(in, out, flags) << std::endl;
    return 0;
  }
  
//  if (debug_mode) {
//    auto false_vacuum = t[jj].false_vacuum;
//    auto true_vacuum = t[jj].true_vacuum;
//    auto TC = t[jj].TC;
//    std::cout << "V^high = " << model.V(false_vacuum, TC) << std::endl;
//    std::cout << "V1T^high = " << model.V1T(false_vacuum, TC) << std::endl;
//    std::cout << "VCW^high = " << model.V1(false_vacuum) << std::endl;
//    std::cout << "daisy^high = " << model.daisy(false_vacuum,TC) << std::endl;
//    std::cout << "Vtree^high = " << model.V0(false_vacuum) << std::endl;
//    
//    std::cout << "V^low = " << model.V(true_vacuum, TC) << std::endl;
//    std::cout << "V1T^low = " << model.V1T(true_vacuum, TC) << std::endl;
//    std::cout << "VCW^low = " << model.V1(true_vacuum) << std::endl;
//    std::cout << "daisy^low = " << model.daisy(true_vacuum,TC) << std::endl;
//    std::cout << "Vtree^low = " << model.V0(true_vacuum) << std::endl;
//    
//    auto d2V_dx2_high = model.d2V_dx2(false_vacuum, TC);
//    std::cout << "d2V_dx2_high = " << d2V_dx2_high << std::endl;
//    auto d2V_dx2_low = model.d2V_dx2(true_vacuum, TC);
//    std::cout << "d2V_dx2_low = " << d2V_dx2_low << std::endl;
//    
//    auto d2V_dxdt_high = model.d2V_dxdt(false_vacuum, TC);
//    std::cout << "d2V_dxdt_high = " << d2V_dxdt_high << std::endl;
//    auto d2V_dxdt_low = model.d2V_dxdt(true_vacuum, TC);
//    std::cout << "d2V_dxdt_low = " << d2V_dxdt_low << std::endl;
//    
//    std::cout << "@ EWSB VEV" << std::endl;
//    Eigen::VectorXd test(2);
//    test <<  SM_parameters[1], 0;
//    std::cout << "dV1T_dT(T=0) = " << model.dV1T_dT(test, 0.0001) << std::endl;
//    std::cout << "ddaisy_dT(T=0) = " << model.ddaisy_dT(test, 0.0001) << std::endl;
//    
//    std::cout << "dV1T_dT(T=TC) = " << model.dV1T_dT(test, TC) << std::endl;
//    std::cout << "ddaisy_dT(T=TC) = " << model.ddaisy_dT(test, TC) << std::endl;
//    
//  }
  
  
  std::vector<double> out = {(float)t.size(), t[jj].TC, t[jj].true_vacuum[0], t[jj].true_vacuum[1], t[jj].false_vacuum[0], t[jj].false_vacuum[1]};
  
  std::cout <<  std::fixed << std::setprecision(4) << model.get_muh_sq() << " & "<< model.get_mus_sq() << " & " << model.get_lambda_h()
  << " & " << t[jj].TC << " & " << t[jj].true_vacuum[0]/t[jj].TC << " & " << t[jj].false_vacuum[1] << " & " <<  t[jj].true_vacuum[0] << std::endl;

  std::cout <<  std::fixed << std::setprecision(8) << model.get_muh_sq() << " & "<< model.get_mus_sq() << " & " << model.get_lambda_h()
  << " & " << t[jj].TC<< " & " << t[jj].true_vacuum[0]/t[jj].TC << " & " << t[jj].false_vacuum[1] << " & " <<  t[jj].true_vacuum[0] << std::endl;
  
  output_file << toString(in, out, flags) << std::endl;
  output_file.close();  
  // Print the data in a particular format for plotting
//  if (debug_mode) PhaseTracer::phase_plotter(tf, "xSM_MSbar");
  return 0;
}
