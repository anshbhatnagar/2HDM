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
#include "omp.h"
#include <nlohmann/json.hpp>
#include "../ProgressBar/progress.hpp"

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


class THDM{
  private:
    json modelParams;
    double mh, mH, mA, mHpm, mu, alpha, beta, rgLam, cosbma, tanb;
    double high_T, low_T;
    double mZ = 91.1876;

  public:
    double trans_alpha;
    double trans_beta_H;
    double trans_Tnuc;

    THDM(json inModelParams){
      modelParams = inModelParams;
    
      LOGGER(fatal);
    
      // Initialise Model
    
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
        exit(0);
      }
    }

    std::string CalcPhaseTransition(){
      std::vector<double> in = {mh, mH, mA, mHpm, mu, alpha, beta, rgLam};
      EffectivePotential::DR_2HDM model(in);

      //model.printLagrangianParams();
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
        return "Couldn't find phases.";
      }
      //std::cout << pf;
    
      // Check there's multiple phases
      auto p1 = pf.get_phases();
      if ( p1.size() == 0 ) { 
        return "Couldn't find multiple phases.";
        }
    
      // Make ActionCalculator object
      PhaseTracer::ActionCalculator ac(model);
      ac.set_action_calculator(PhaseTracer::ActionMethod::BubbleProfiler);
      
      // Make TransitionFinder object and find the transitions
      PhaseTracer::TransitionFinder tf(pf, ac);
      tf.find_transitions();
    
      //std::cout << tf;
    
      auto t = tf.get_transitions();
      if (t.size()==0){
        return "Couldn't find transitions.";
      }
    
      //PhaseTracer::potential_plotter(model,t[0].TC,"2HDM");
      //PhaseTracer::potential_line_plotter(model,t,"2HDM");
    
      if ( isnan(t[0].TN) ){
        return "Couldn't find finite nucleation temperature.";
      }
    
      trans_Tnuc = t[0].TN;
    
      // Make GravWave Object
      PhaseTracer::GravWaveCalculator gc(tf);
      gc.set_min_frequency(1e-4);
      gc.set_max_frequency(1e+1);
      gc.calc_spectrums();
    
      // std::cout << gc;
    
      auto gw = gc.get_spectrums();
    
      if ( isnan(gw[0].alpha) || isnan(gw[0].beta_H)){
        return "Couldn't find gravitational waves.";
      }
      trans_alpha = gw[0].alpha;
      trans_beta_H = gw[0].beta_H;

      return "Success!";
      
    }
};



int main(int argc, char* argv[]){
  json runParams;
  std::vector<json> modelParamsVec;
  std::string runVar;
  int numRuns;
  json output;
  std::string outputFileName;
  
  output = json::array();

  try{
    runParams = readFile("runParams.json");
    runVar = runParams["Run Variable"].get<std::string>();
    numRuns = runParams["Number of Runs"].get<int>();
    outputFileName = runParams["Output File Name"].get<std::string>();

    double startVar = runParams[runVar]["Start"].get<double>();
    double endVar = runParams[runVar]["End"].get<double>();
    assert(endVar > startVar);
    double increVar = (endVar - startVar)/numRuns;

    for(int i = 0; i < numRuns; i++){
      
      json modelParams = runParams;
      modelParams.erase(modelParams.find("Run Variable"));
      modelParams.erase(modelParams.find("Number of Runs"));
      modelParams.erase(modelParams.find(runVar));
      modelParams[runVar] = startVar + i*increVar;
      

      modelParamsVec.push_back(modelParams);
    }
    
  } catch(...){
    std::cout << "Input file might not be formatted correctly! Quitting now." << std::endl;
    return 0;
  }

  progress::progBar bar(5*numRuns);

  std::thread progress;
  progress = bar.start();

  #pragma omp parallel for shared(modelParamsVec, output, bar)
    for(int i = 0; i < numRuns; i++){
      json runOut;
      std::string runError;
      int errCount = 0;
      int numTries = 5;
      json testParams = modelParamsVec[i];

      runOut[runVar] = modelParamsVec[i][runVar].get<double>();

      //omp_set_num_threads(4);

      while(errCount < numTries){

        try{
        THDM thdm(testParams);
        
        runError = thdm.CalcPhaseTransition();

        if(runError == "Success!"){
          runOut["alpha"] = thdm.trans_alpha;
          runOut["beta/H"] = thdm.trans_beta_H;
          runOut["Tnuc"] = thdm.trans_Tnuc;
        } else if(runError == "Couldn't find gravitational waves."){
          runOut["Tnuc"] = thdm.trans_Tnuc;
          runOut["Error"] = runError;
        } else {
          runOut["Error"] = runError;
        }
  
        #pragma omp critical
        {
          output.push_back(runOut);
          bar.add(numTries - errCount);
        }
        
        errCount = numTries;

        } catch(...){
          if(errCount == numTries - 1){
            runOut["Error"] = "Failed to run a parameter point!";
            #pragma omp critical
            {
              output.push_back(runOut);
              bar.add(numTries - errCount);
            }
            errCount = numTries;
          } else{
          testParams[runVar] = 1.001*testParams[runVar].get<double>();
          runOut[runVar] = testParams[runVar].get<double>();
          errCount += 1;
          #pragma omp critical
          {
            bar.add(1);
          }
        }
        }
     }
    
  }
    //std::cout<<output << std::endl;

    bar.finish();
    progress.join();

    std::ofstream outFile;
    outFile.open(outputFileName+".json");

    outFile << std::setw(4) << output << std::endl;

    outFile.close();

    std::cout<<std::endl;

    return 0;

    


}