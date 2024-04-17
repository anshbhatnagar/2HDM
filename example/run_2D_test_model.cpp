/**
  2D example program for PhaseTracer.
*/

#include <iostream>

#include "models/2D_test_model.hpp" // Located in effective-potential/include/models
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "logger.hpp"


int main(int argc, char* argv[]) {

  const bool debug_mode = argc > 1 && strcmp(argv[1], "-d") == 0;

  // Set level of screen  output
  if (debug_mode) {
      LOGGER(debug);
  } else {
      LOGGER(fatal);
  }

  // Construct our model
  EffectivePotential::TwoDimModel model;

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.find_phases();
  std::cout << pf;
  
  // Make PhaseFinder object
  PhaseTracer::ActionCalculator ac(model);
  ac.set_use_BubbleProfiler(false);
  
  double T=77.61682;
//  auto vacuums = pf.find_minima_at_t(T);
  
  Eigen::VectorXd true_vacuum(2);
  true_vacuum << 289.27693787, 389.92593071;
  Eigen::VectorXd false_vacuum(2);
  false_vacuum <<  234.2959619,  -111.48232561;
  
  ac.get_action(true_vacuum,false_vacuum,T);
  
//  // Make TransitionFinder object and find the transitions
//  PhaseTracer::TransitionFinder tf(pf,ac);
//  tf.find_transitions();
//  std::cout << std::setprecision (15) << tf;

//  if (debug_mode) {
//    PhaseTracer::potential_plotter(model, tf.get_transitions().front().TC, "2D_test_model", 0., 2., 0.01, -2., 0., 0.01);
//    PhaseTracer::potential_line_plotter(model, tf.get_transitions(), "2D_test_model");
//    PhaseTracer::phase_plotter(tf, "2D_test_model");
//  }
}
