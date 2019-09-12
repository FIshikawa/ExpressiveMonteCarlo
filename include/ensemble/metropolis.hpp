#ifndef ENSEMBLE_METROPOLIS_HASTING_HPP
#define ENSEMBLE_METROPOLIS_HASTING_HPP

#include <utility>
#include <cmath>
#include <string>
#include <random>
#include <functional>

namespace ensemble{

template<typename Potential, typename Proposer>
class Metropolis{
public:
  static std::string name() { return "Metropolis method"; }
  Metropolis(Potential target_potential,Proposer proposer) 
            {target_potential_ = target_potential; proposer_ = proposer;}
  void montecalro(std::vector<double>& z, std::mt19937 & mt, double temperture, int& counter){
    // const double kB = 1.38064852 / pow(10.0,23.0);
    std::pair<int, double> new_values = proposer_(z,mt);
    int target = new_values.first;
    double past_value = z[target];
    double new_value =  new_values.second;
    z[target] = new_value;
    double ene_new = target_potential_(target,z);
    z[target] = past_value;
    double ene_past = target_potential_(target,z);
    
    double acceptance = std::exp(-1.0/temperture * ene_new)/std::exp(-1.0/temperture * ene_past);
    std::uniform_real_distribution<> uniform_random(0,1.0);
    double dice = uniform_random(mt);
   
    if(acceptance > dice){
      counter += 1;
      z[target] = new_value;
    }
  }

private:
  Potential target_potential_;
  Proposer proposer_; 
};

} //end namespace

#endif //ENSEMBLE_METROPOLIS_HPP

