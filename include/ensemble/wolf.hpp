#ifndef ENSEMBLE_WOLF_HPP
#define ENSEMBLE_WOLF_HPP

#include <cmath>
#include <string>
#include <random>
#include <stack>

namespace ensemble{

class Wolf{
public:
  static std::string name() { return "Single cluster update, Wolf algorithm"; }
  Wolf(double J, int num) : J_(J), num_(num) {}

  template <class Hamiltonian, class Proposer>
  void montecalro(std::vector<double>z, Hamiltonian const& hamiltonian, Proposer const& proposer, std::mt19937 & mt, double temperture, int& counter){
    std::uniform_int_distribution<> dist(0,num_-1);
    std::uniform_real_distribution<> Spin_Rand(0,2*M_PI);
    std::uniform_real_distribution<> realdist(0,1.0);
    std::stack<int> spin_checklist;
    double *spin = &z[0];
    int check_tag = dist(mt);
    double phi = Spin_Rand(mt); //project angle
    int spin_original = sign(std::cos(spin[check_tag] - phi));
    spin[check_tag] = -spin[check_tag] + 2.*phi + M_PI;
    spin_checklist.push(check_tag);
    int cluster_size = 0;
    while(!spin_checklist.empty()){
      ++cluster_size;
      int center_tag = spin_checklist.top();
      double spin_center = std::cos(spin[center_tag] - phi);
      spin_checklist.pop();
      int Nd_t = hamiltonian.Nd();
      for(int k = 0; k < Nd_t; ++k){
        int nearest_tag = hamiltonian.table(center_tag,k);
        double spin_nearest = std::cos(spin[nearest_tag] - phi);
        double J_int = J_ * std::fabs(spin_nearest * spin_center);
        spin_nearest = sign(spin_nearest);
        double P_accept = 1 - std::exp(-2.*J_int/temperture);
        if(spin_nearest == spin_original && realdist(mt) < P_accept){
          spin_checklist.push(nearest_tag);
          spin[nearest_tag] = -spin[nearest_tag] + 2.*phi + M_PI;
        }
      }
    } 
    counter += cluster_size;
  }

private:
  int sign(double hoge)const{ return (hoge > 0) - (hoge < 0);} 
  int num_;
  double J_;
}; //equil end

} //end namespace

#endif //ENSEMBLE_WOLF_CLASSICAL_XY_HPP
