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

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "thread"
#include <nlohmann/json.hpp>

#include "DRalgo_2HDM.hpp"
#include "phasetracer.hpp"

using json = nlohmann::json;

json readFile(std::string fileName){
  std::ifstream file;
  file.open(fileName);

  if(file.fail()){
      throw std::runtime_error("error loading model parameters file!");
  }

  json data = json::parse(file);

  file.close();

  return data;
}

std::string toString(std::vector<double> in, std::vector<double> out) {
  std::stringstream data_str;
  for (auto i : in    ) data_str << i << "\t";
  for (auto i : out   ) data_str << i << "\t";
  return data_str.str();
}


int main(int argc, char* argv[]) {

  std::ofstream outFile;
  outFile.open("output.json");

  LOGGER(debug);

  double mh, mH, mA, mHpm, mu, alpha, beta, rgLam, cosbma, tanb;
  double high_T, low_T;
  double mZ = 91.1876;

  // Initialise Model

  json modelParams = readFile("modelParams.json");

  try {
    cosbma = modelParams["cos(beta - alpha)"].get<double>();
    tanb = modelParams["tan(beta)"].get<double>();

    mh = modelParams["mh"].get<double>();
    mH = modelParams["mH"].get<double>();
    mA = modelParams["mA"].get<double>();
    mHpm = modelParams["mHpm"].get<double>();
    mu = modelParams["mu"].get<double>();;
    beta = atan(tanb);
    alpha = beta - acos(cosbma);
    rgLam = modelParams["Input Scale"].get<double>();

    high_T = modelParams["High T"].get<double>();
    low_T = modelParams["Low T"].get<double>();
  } catch(...){
    std::cout << "Input file might not be formatted correctly! Quitting now." << std::endl;
    return 0;
  }

  std::vector<double> in = {mh, mH, mA, mHpm, mu, alpha, beta, rgLam};
  EffectivePotential::DR_2HDM model(in);

  model.printLagrangianParams();
  //model.printHiggsvev();

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  
  pf.set_seed(0);
  pf.set_check_hessian_singular(false);
  pf.set_check_vacuum_at_high(false);
  pf.set_t_high(high_T);
  pf.set_t_low(low_T);

  try {
    pf.find_phases();
  } catch (...) {
    std::vector<double> out = {-1, 0, 0, 0, 0, 0, 0, 0};
    outFile << toString(in, out) << std::endl;
    return 0;
  }
  std::cout << pf;

  // Check there's multiple phases
  auto p1 = pf.get_phases();
  if ( p1.size() == 0 ) { return 0; }

  // Make ActionCalculator object
  PhaseTracer::ActionCalculator ac(model);
  //ac.set_action_calculator(PhaseTracer::ActionMethod::BubbleProfiler);
  
  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf, ac);

  json output = tf.get_S_T(100);

  outFile << std::setw(4) << output << std::endl;

  outFile.close();

  tf.find_transitions();

  std::cout << tf;

  std::cout<<std::endl;

  return 0;
  
}