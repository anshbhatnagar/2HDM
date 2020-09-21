/**
  Example of h-bar expansion using xSM model.
*/

#include <iostream>

#include "models/xSM_MSbar.hpp"
#include "transition_finder.hpp"
#include "h_bar_expansion.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "logger.hpp"


int main(int argc, char* argv[]) {
  LOGGER(debug);

  const double lambda_hs = atof(argv[1]);
  const double Q = atof(argv[2]);
  std::cout << "lambda_hs = " << lambda_hs << std::endl
            << "Q = " <<  Q << std::endl;

  // Construct our model
  auto model = EffectivePotential::make_xSM(lambda_hs, Q);
  model.set_daisy_method(EffectivePotential::DaisyMethod::None);

  // Make PhaseFinder object and find the phases
  PhaseTracer::HbarExpansion hb(model);
  hb.set_seed(0);
  hb.find_phases();
  std::cout << hb;

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(hb);
  tf.find_transitions();
  std::cout << std::setprecision(15) << tf;
  return 0;
}
