#ifndef ENSEMBLE_HYBRID_MONTECARLO_HPP
#define ENSEMBLE_HYBRID_MONTECARLO_HPP

#include <cmath>
#include <string>
#include <random>
#include <functional>

namespace ensemble{

template<typename Potential, class Integrator, typename Gradient>
class HybridMonteCarlo{
public:
  static std::string name() { return "Hybrid Monte Calro method"; }
  HybridMonteCarlo(Potential potential, Gradient gradient, int num_particles, int time_steps, double dt)
                  : num_particles_(num_particles),integrator_(2*num_particles),k1_(2*num_particles),k2_(2*num_particles), time_steps_(time_steps),dt_(dt)
                  {potential_ = potential; force_set(potential,gradient);}

  void montecalro(std::vector<double>& z, std::mt19937 & mt, double temperture, int& counter){
    std::normal_distribution<> normal_dist(0.0,std::sqrt(temperture));
    for(int i = 0;i < num_particles_;++i) k1_[i] = z[i];
    for(int i = 0;i < num_particles_;++i) k1_[i+num_particles_] = normal_dist(mt);
    k2_ = k1_;
    for(int step = 0 ; step < time_steps_;++step) integrator_.step(0.0, dt_, k2_, force_);
    double ene_past = potential_(k1_) + kinetic(k1_);
    double ene_new  = potential_(k2_) + kinetic(k2_);
    double acceptance = std::exp(-1.0/temperture * ene_new)/std::exp(-1.0/temperture * ene_past);
    std::uniform_real_distribution<> uniform_random(0,1.0);
    double dice = uniform_random(mt);
    if(acceptance > dice){
      counter += 1;
      for(int i = 0;i < num_particles_;++i) z[i] = k2_[i];
    }
  }

  double kinetic(std::vector<double> z){
    double ene = potential_(z);
    for(int i = 0;i < num_particles_;++i) ene += 0.5 * z[num_particles_+i] * z[num_particles_+i];
    return ene;
  }

  void force_set(Potential potential, Gradient gradient){
    int num_particles = num_particles_;
    force_ = [num_particles, gradient, potential](double t, std::vector<double> & z, std::vector<double>& force){
      double *x = &z[0];
      double *v = &z[num_particles];
      double *fx = &force[0];
      double *fv = &force[num_particles];
      for(int i = 0; i < num_particles; ++i) fx[i] = v[i];
      for(int i = 0; i < num_particles ; ++i) fv[i] = gradient(i,z,potential);
    };
  }

private:
  Potential potential_;
  Integrator integrator_; 
  std::function<void(double,std::vector<double>&,std::vector<double>&)> force_;
  int num_particles_,time_steps_;
  double dt_;
  mutable std::vector<double> k1_,k2_;
};

} //end namespace

#endif //ENSEMBLE_HYBRID_MONTECARLO_HPP
