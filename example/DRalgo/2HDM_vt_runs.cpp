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
    double v_T;

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
        mu = pow(10.,modelParams["log10(mu)"].get<double>());
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
      //ac.set_action_calculator(PhaseTracer::ActionMethod::BubbleProfiler);

      // Make TransitionFinder object and find the transitions
      PhaseTracer::TransitionFinder tf(pf, ac);


      //std::cout << tf;
      try{
        v_T = tf.get_v_T();
      }
      catch(...){
        return "Failed to get v/T!";
      }

    
      
      return "Success!";
      
    }
};



int main(int argc, char* argv[]){
  json runParams;
  std::vector<std::vector<json>> modelParamsVec;
  std::string runVarX;
  std::string runVarY;
  int numRunsX = 1;
  int numRunsY = 1;
  json output;
  std::string outputFileName;

  bool twoDim = false;
  
  output = json::array();

  try{
    runParams = readFile("runParams.json");

    double startVarX;
    double startVarY;
    double endVarX;
    double endVarY;

    if(runParams.contains("Run Variable")){
      runVarX = runParams["Run Variable"].get<std::string>();
      numRunsX = runParams["Number of Runs"].get<int>();
      startVarX = runParams[runVarX]["Start"].get<double>();
      endVarX = runParams[runVarX]["End"].get<double>();
      
      runParams.erase(runParams.find("Run Variable"));
      runParams.erase(runParams.find("Number of Runs"));
      runParams.erase(runParams.find(runVarX));
    }else if(runParams.contains("Run Variable Y")){
      twoDim = true;
      runVarX = runParams["Run Variable X"].get<std::string>();
      runVarY = runParams["Run Variable Y"].get<std::string>();
      numRunsX = runParams["Number of Runs in X"].get<int>();
      numRunsY = runParams["Number of Runs in Y"].get<int>();
      startVarX = runParams[runVarX]["Start"].get<double>();
      endVarX = runParams[runVarX]["End"].get<double>();
      startVarY = runParams[runVarY]["Start"].get<double>();
      endVarY = runParams[runVarY]["End"].get<double>();

      runParams.erase(runParams.find("Run Variable X"));
      runParams.erase(runParams.find("Number of Runs in X"));
      runParams.erase(runParams.find(runVarX));
      runParams.erase(runParams.find("Run Variable Y"));
      runParams.erase(runParams.find("Number of Runs in Y"));
      runParams.erase(runParams.find(runVarY));
    }else{
      throw 0;
    }
    outputFileName = runParams["Output File Name"].get<std::string>();

    
    assert(endVarX > startVarX);
    double increVarX = (endVarX - startVarX)/numRunsX;
    double increVarY;
    if(twoDim){
      assert(endVarY > startVarY);
      increVarY = (endVarY - startVarY)/numRunsY;
    }


    if(!twoDim){
      std::vector<json> modelParams1D;
      for(int i = 0; i < numRunsX; i++){
        json modelParams = runParams;
        modelParams[runVarX] = startVarX + i*increVarX;
  
        modelParams1D.push_back(modelParams);
      }
      modelParamsVec.push_back(modelParams1D);

    }else{
      for(int i = 0; i < numRunsX; i++){
        std::vector<json> modelParams1D;
        for(int j = 0; j < numRunsY; j++){
          json modelParams = runParams;

          modelParams[runVarX] = startVarX + i*increVarX;

          modelParams[runVarY] = startVarY + j*increVarY;
    
          modelParams1D.push_back(modelParams);
        }
        modelParamsVec.push_back(modelParams1D);
      }
    }
    
  } catch(...){
    std::cout << "Input file might not be formatted correctly! Quitting now." << std::endl;
    return 0;
  }
  
  int numTries = 5;

  progress::progBar bar(numTries*numRunsX*numRunsY);

  std::thread progress;
  progress = bar.start();

  #pragma omp parallel for shared(modelParamsVec, output, bar)
  for(int i = 0; i < numRunsX; i++){
    for(int j = 0; j < numRunsY; j++){
      json runOut;
      std::string runError;
      int errCount = 0;
      json testParams = modelParamsVec[i][j];

      
      runOut[runVarX] = modelParamsVec[i][j][runVarX].get<double>();

      if(twoDim){
        runOut[runVarY] = modelParamsVec[i][j][runVarY].get<double>();
      }

      //omp_set_num_threads(4);

      while(errCount < numTries){

        try{
        THDM thdm(testParams);
        
        runError = thdm.CalcPhaseTransition();
        
        if(runError != "Success!"){
          runOut["Error"] = runError;
        }else{
          runOut["v_T"] = thdm.v_T;
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
            testParams[runVarX] = 1.0001*testParams[runVarX].get<double>();
            runOut[runVarX] = testParams[runVarX].get<double>();
            if(twoDim){
              testParams[runVarY] = 1.0001*testParams[runVarY].get<double>();
              runOut[runVarY] = testParams[runVarY].get<double>();
            }
            
            errCount += 1;
            #pragma omp critical
            {
              bar.add(1);
            }
          }
        }
     }
    
    }
  }

  bar.finish();
  progress.join();

  std::ofstream outFile;
  outFile.open(outputFileName+".json");

  outFile << std::setw(4) << output << std::endl;

  outFile.close();

  std::cout<<std::endl;

  return 0;

    


}