#include <gtest/gtest.h>
#include <physics/harmonic_oscillator_fixed_end.hpp>
#include <lattice/chain_open_boundary.hpp>
#include <integrator/runge_kutta_4th.hpp>

namespace {
class HarmonicOscillatorFixedEndTest: public ::testing::Test {
  protected:  
  virtual void SetUp(){
    // set by lattice::Chain
    int Ns = 10;
    lattice::ChainOpenBoundary lattice(Ns);
    num_particles = lattice.set_num_particles(Ns); 
    ASSERT_EQ(num_particles, Ns);
    N_adj = lattice.number_adjacent();
    ASSERT_EQ(N_adj,2);
    pair_table = std::vector<std::vector<int> >(num_particles,std::vector<int>(N_adj));
    lattice.create_table(pair_table);
    ASSERT_EQ(pair_table[0][0],0);
    ASSERT_EQ(pair_table[0][1],1);
    ASSERT_EQ(pair_table[num_particles-1][0],num_particles - 2);
    ASSERT_EQ(pair_table[num_particles-1][1],num_particles - 1);
    // set interaction constants 
    J = 2.0;
    // set hamiltonian
    z.resize(2*num_particles);
    for(int i = 0 ; i < num_particles ; ++i){
      z[i+num_particles] = 1.0;
      if(i % 2 == 0) z[i] = 1.0;
      else z[i] = 0;
    }
    expected_energy_kinetic = 0.5 * num_particles;
    expected_energy_potential = 1.0 * num_particles;
    expected_energy_total = expected_energy_kinetic + expected_energy_potential;
  }

  int  num_particles, N_adj;
  double J;
  double expected_energy_potential, expected_energy_kinetic, expected_energy_total;
  std::vector<std::vector<int> > pair_table;
  std::vector<double> z;
};

TEST_F(HarmonicOscillatorFixedEndTest, EnergyTest) {
  hamiltonian::HarmonicOscillatorFixedEnd hamiltonian(num_particles,J,pair_table,N_adj);
  double energy_potential = 0.0;
  double energy_kinetic = 0.0;
  for(int i = 0; i < num_particles; ++i){
    energy_potential += hamiltonian.target_potential_energy(i,z,0);
    energy_kinetic += hamiltonian.target_kinetic_energy(i,z);
  }
  energy_potential += hamiltonian.target_potential_energy(-1,z,0);
  energy_potential *= 0.5;
  EXPECT_DOUBLE_EQ(energy_kinetic, hamiltonian.kinetic_energy(0,z));
  EXPECT_DOUBLE_EQ(energy_potential, hamiltonian.potential_energy(0,z));
  EXPECT_DOUBLE_EQ(expected_energy_potential, hamiltonian.potential_energy(0,z));
  EXPECT_DOUBLE_EQ(expected_energy_kinetic, hamiltonian.kinetic_energy(0,z));
  EXPECT_DOUBLE_EQ(expected_energy_total, hamiltonian.energy(0,z));

}

TEST_F(HarmonicOscillatorFixedEndTest, ForceTest) {
  hamiltonian::HarmonicOscillatorFixedEnd hamiltonian(num_particles,J,pair_table,N_adj);
  // set integrator 
  integrator::RungeKutta4th integrator(2*num_particles);
  double dt = 1.0e-2;
  double pt = 0.0;
  double N_time = 1e+4;
  double error_rungekutta = dt * dt * dt * dt * dt * 10; //10 is assumed constant of error 
  double energy_init = hamiltonian.energy(pt,z);
  double potential_init = hamiltonian.potential_energy(pt,z);
  double kinetic_init = hamiltonian.kinetic_energy(pt,z);
  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z,hamiltonian);
    pt += dt;
  }
  EXPECT_FALSE(energy_init == hamiltonian.energy(pt,z));
  EXPECT_NEAR(energy_init/hamiltonian.energy(pt,z),1.0,error_rungekutta*N_time);
  EXPECT_FALSE(potential_init == hamiltonian.potential_energy(pt,z)); 
  EXPECT_FALSE(kinetic_init == hamiltonian.kinetic_energy(pt,z)); 
  for(int i = 0; i < N_time; ++i){
    integrator.step(pt,dt,z,hamiltonian);
    pt += dt;
  }
  EXPECT_FALSE(energy_init == hamiltonian.energy(pt,z));
  EXPECT_NEAR(energy_init/hamiltonian.energy(pt,z),1.0,error_rungekutta*N_time);
  EXPECT_FALSE(potential_init == hamiltonian.potential_energy(pt,z)); 
  EXPECT_FALSE(kinetic_init == hamiltonian.kinetic_energy(pt,z)); 
}

}//end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}