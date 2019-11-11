#ifndef HAMILTONIAN_HARMONICALLY_TRAPPED_ATOM_HPP
#define HAMILTONIAN_HARMONICALLY_TRAPPED_ATOM_HPP

#include <string>
#include <cmath>
#include <vector>

namespace hamiltonian{

class HarmonicallyTrappedAtom{
public:
  static std::string name() { return "Harmonically trapped atom"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  // H = 0.5*(v^2 + w^2 + ...) + 0.5*J*(x^2 + y^2 + ...)
  HarmonicallyTrappedAtom(
                          int num_particles, 
                          double J, 
                          std::vector<std::vector<int > > table, 
                          int num_adjacent
                          ) : num_particles_(num_particles), 
                              J_(J), 
                              table_(table),
                              nd_(num_adjacent){}

  double energy(const double t, std::vector<double> const& z) const {
    return potential_energy(t, z) + kinetic_energy(t, z);
  }
  double potential_energy(const double  t, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_ ; ++i) ene += target_potential_energy(i,z,t);
    return ene/2.0;
  }
  double kinetic_energy(const double t, std::vector<double> const& z) const {
    double ene = 0;
    for(int i = 0; i < n_; ++i) ene += target_kinetic_energy(i,z);
    return ene;
  }

  double target_potential_energy(int l, std::vector<double> const& z, double t) const {
    const double *x = &z[0];
    return  0.5 * J_ * x[l] * x[l];
  }

  double target_kinetic_energy(int l, std::vector<double> const& z) const {
    const double *v = &z[n_];
    return 0.5 * v[l] * v[l];
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
        fv[i] -= J_ * x[i];
      }
    }
  }

  int nd()const { return nd_ ;}
  int table(int i, int j)const{ return table_[i][j];} 

private:
  int n_;
  double j_;
  std::vector<std::vector<int> > table_;
  int nd_;
};

} //end namespace

#endif //HAMILTONIAN_HARMONICALLY_TRAPPED_ATOM_HPP

