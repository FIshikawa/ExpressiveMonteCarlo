#ifndef INTEGRATOR_VELOCITY_VELRET_HPP
#define INTEGRATOR_VELOCITY_VELRET_HPP

#include <string>
#include <vector>

namespace integrator {

class VelocityVelret{
public:
  static std::string name() { return "Velocity Verlet method"; }
  // z[0...n-1]: position, z[n...2n-1]: momentum
  VelocityVelret(unsigned int dim) : dim_(dim), k1_(dim), k2_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    const double h22 = h * h / 2;
    const double h2 = h / 2;
    const int  dim2 = dim_ / 2;
    f(t, y, k1_);
    for (int i = 0; i < dim2; ++i) y[i] = y[i] +  h * k1_[i] + h22 * k1_[i + dim2];
    f(t + h, y, k2_);
    for (int i = dim2; i < dim_; ++i) y[i] = y[i] +  h2 * (k2_[i]+k1_[i]);
  }
private:
  unsigned int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
};

} // namespace integrator

#endif // INTEGRATOR_VELOCITY_VELRET_HPP
