#ifndef HAMILTONIAN_HARMONIC_OSCILLATOR_FIXED_END_HPP
#define HAMILTONIAN_HARMONIC_OSCILLATOR_FIXED_END_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class HarmonicOscillatorFixedEnd  {
public:
  static std::string name() { return "Classical Harmonic Oscillator : Fixed End : must use Open Boundary Chain Lattice"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  HarmonicOscillatorFixedEnd(int num_particles, double J, std::vector<std::vector<int > > table, int N_adj) : n_(num_particles), J_(J), table_(table),Nd_(N_adj){}
  double energy(double t, std::vector<double> const& z) const {
    return potential_energy(t, z) + kinetic_energy(t, z);
  }
  double potential_energy(double  t, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = -1; i < n_ ; ++i) ene += target_potential_energy(i,z,t);
    return ene/2.0;
  }
  double kinetic_energy(double /* t */, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_; ++i) ene += target_kinetic_energy(i,z);
    return ene;
  }

  double target_potential_energy(int l, std::vector<double> const& z, double t) const {
    const double *x = &z[0];
    double ene = 0;
    if( l != -1){
      if(l == 0 || l == n_ - 1) ene += 0.5 * J_ * x[l] * x[l];
      for(int i = 0; i < Nd_ ; ++i){
        double dx = x[l] - x[table_[l][i]];
        ene += 0.5* J_ * dx * dx;
      }
    }
    else for(auto& i : {0, n_ -1}) ene += 0.5 * J_ * x[i] * x[i];
    return ene;
  }

  double target_kinetic_energy(int l, std::vector<double> const& z) const {
    const double *v = &z[n_];
    double ene = 0;
    ene += 0.5 * v[l]*v[l];
    return ene;
  }

  // "force" calculation
  void operator()(double t, std::vector<double> const& z, std::vector<double>& force) const {
    const double *x = &z[0];
    const double *v = &z[n_];
    double *fx = &force[0];
    double *fv = &force[n_];
    for(int i = 0; i < n_ ; ++i) fx[i] = v[i];
    for(int i = 0; i < n_ ; ++i){
      fv[i] = 0.0;
      for(int d = 0; d < Nd_; ++d){
        fv[i] -=  J_ * (x[i] -x[table_[i][d]]);
      }
    }
    for(auto& i : {0, n_ - 1}) fv[i] -= J_ * x[i];
  }

  int Nd()const { return Nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

private:
  int n_;
  double J_;
  std::vector<std::vector<int> > table_;
  int Nd_;
};


} //end namespace

#endif //HAMILTONIAN_HARMONIC_OSCILLATOR_FIXED_END_HPP 
