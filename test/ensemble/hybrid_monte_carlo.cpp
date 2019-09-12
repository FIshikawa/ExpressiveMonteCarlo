#include <gtest/gtest.h>
#include <limits>
#include <lattice/chain.hpp>
#include <ensemble/hybrid_monte_carlo.hpp>
#include <integrator/position_velret.hpp>

using Potential = std::function<double(std::vector<double> const &)>;
using Gradient = std::function<double(int,std::vector<double>&,Potential)>; 
using Integrator = integrator::PositionVelret;

namespace {
class HybridMonteCarloEnsemblerTest : public ::testing::Test {
protected:
  virtual void SetUp(){

    // set constants
    num_particles = 1;
    J = 1.0;
    precision = std::numeric_limits<double>::epsilon();
    temperture = 1.0;
    time_steps = 1e+2;
    dt = 1e-2;

    // set vector
    z.resize(num_particles);
    for(int i = 0; i < num_particles; ++i) z[i] = 0.0;

    // set functions
    // set potential
    double J_temp = J;
    int num_particles_temp = num_particles;
    potential = [J_temp,num_particles_temp](std::vector<double> const & z){
      double ene = 0;
      for(int i = 0; i < num_particles_temp; ++i) ene += J_temp * z[i] * z[i]; 
      return ene;
    };

    // set gradient
    double epsilon = 1e-3;
    gradient = [epsilon](int i, std::vector<double> & z, Potential potential){
      double *x = &z[0];
      double x_temp = x[i];
      double x_backward = x[i] - epsilon;
      double x_forward = x[i] + epsilon;
      x[i] = x_backward;
      double ene_backward = potential(z);
      x[i] = x_forward;
      double ene_forward= potential(z);
      x[i] = x_temp;
      return -1.0*(ene_forward - ene_backward)/(2.0*epsilon);
    };
    // alternartively, 
    // gradient = [J_temp](int i , std::vector<double> & z, Potential potential){
    //   double *x = &z[0];
    //   return -1.0 * J_temp * x[i]
    // }
  } 
  
  Potential potential;
  Gradient gradient;
  std::vector<double> z;
  int num_particles,time_steps;
  double J,precision,temperture,dt;
  std::vector<std::vector<int> > pair_table;
};

TEST_F(HybridMonteCarloEnsemblerTest, BasicTest){
  int counter = 0;
  std::size_t seed = 1234;
  std::mt19937 mt(seed);

  Integrator integrator(2*num_particles);
  ensemble::HybridMonteCarlo<Potential, Integrator, Gradient> hybrid_monte_carlo(potential,gradient,num_particles,time_steps,dt);
  hybrid_monte_carlo.montecalro(z,mt,temperture,counter);
  ASSERT_EQ(counter,1);
} 

}// end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
