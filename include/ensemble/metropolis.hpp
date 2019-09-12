#ifndef ENSEMBLE_METROPOLIS_HASTING_HPP
#define ENSEMBLE_METROPOLIS_HASTING_HPP

#include <utility>
#include <cmath>
#include <string>
#include <random>
#include <functional>

namespace ensemble{

class Metropolis{
public:
  static std::string name() { return "Metropolis method"; }
  Metropolis(
             std::function<double(int,std::vector<double>,double)> target_potential,
             std::function<std::pair<int,double>(std::vector<double>,std::mt19937)> proposer
             ) : target_potential_(target_potential),proposer_(proposer){}
  void montecalro(std::vector<double>z, std::mt19937 & mt, double temperture, int& counter){
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::pair<int, double> new_values = proposer_(z,mt);
    int target = new_values.first;
    double dP = 0.0;
    double P_accept;
    double past_value = z[target];
    double new_value =  new_values.second;
    z[target] = new_value;
    double new_dE = target_potential_(target,z,0.0);
    z[target] = past_value;
    double past_dE = target_potential_(target,z,0.0);
    
    dP = std::exp(-1.0/temperture * new_dE)/std::exp(-1.0/temperture * past_dE);
    std::uniform_real_distribution<> uniform_random(0,1.0);
    P_accept = uniform_random(mt);
   
    if(dP > P_accept){
      counter += 1;
      z[target] = new_value;
    }
  }

private:
  std::function<double(int,std::vector<double>,double)> target_potential_;
  std::function<std::pair<int,double>(std::vector<double>,std::mt19937)> proposer_; 
};

} //end namespace

#endif //ENSEMBLE_METROPOLIS_HPP

