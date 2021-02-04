#include <gtest/gtest.h>
#include <limits>
#include <lattice/chain.hpp>
#include <ensemble/metropolis.hpp>


using Potential = std::function<double(int,std::vector<double> const &)>;
using Proposer = std::function<std::pair<int,double>
                          (std::vector<double> const &,std::mt19937 &)>;

namespace {
class MetropolisEnsemblerTest : public ::testing::Test {
protected:
  virtual void SetUp(){
    // set lattice 
    int system_size = 5;
    lattice::Chain lattice(system_size);
    num_particles = lattice.set_num_particles(system_size); 
    num_adjacent = lattice.number_adjacent();
    pair_table = std::vector<std::vector<int> >
                    (num_particles,std::vector<int>(num_adjacent));
    lattice.create_table(pair_table);

    // set constants
    J = -1.0;
    temperture = 3;

    // set vector
    z.resize(num_particles);
    for(int i = 0; i < num_particles; ++i) z[i] = std::pow(1,i);

    // set functions
    // set target potential
    int num_particles_temp = num_particles;
    proposer = 
      [num_particles_temp](std::vector<double> const & z, std::mt19937 & mt){
        std::uniform_int_distribution<> index_dist(0,num_particles_temp-1);
        int target = index_dist(mt);
        double new_value = ((z[target] < 0) ? 1 : -1);
        std::pair<int,double> updates = std::make_pair(target, new_value);
        return updates;
      };

    // set proposer 
    int num_adjacent_temp = num_adjacent;
    double J_temp = J;
    std::vector<std::vector<int> > pair_table_temp = pair_table;
    std::cout << "num_adjacent_temp : " << num_adjacent_temp << std::endl;
    for(int i = 0; i < num_particles; ++i){
      for(int j = 0; j < num_adjacent; ++j)
        std::cout << " i : j -> " << i << " : " << j << " pair : "
                  << pair_table_temp[i][j] << std::endl;
    }
    target_potential = 
      [num_adjacent_temp,J_temp,pair_table_temp]
      (int target, std::vector<double> const & z){
        double ene = 0; 
        for(int i = 0 ; i < num_adjacent_temp; ++i) 
          ene += J_temp*z[target]*z[pair_table_temp[target][i]];
        return ene;
      };
  } 

  Proposer proposer; 
  Potential target_potential;
  std::vector<double> z;
  int num_particles, num_adjacent;
  double J,temperture;
  std::vector<std::vector<int> > pair_table;
};

TEST_F(MetropolisEnsemblerTest, BasicTest){
  std::size_t seed = 42;
  std::mt19937 mt(seed);

  std::vector<double> z_previous(z);
  ensemble::Metropolis<Potential, Proposer> 
                        metropolis(target_potential,proposer);
  int num_iteration = 1;
  for(int counter = 0; counter < num_iteration;){
    std::cout << "<< counter : " << counter << " >>" << std::endl;
    std::cout << "  [before]" << std::endl;
    for(int i = 0; i < num_particles; ++i){
      std::cout << "    i : " << i << " z[i] : " << z[i] << std::endl;
      std::cout << "    target E : " 
                << target_potential(i,z) 
                << std::endl;
    }
    metropolis.montecalro(z,mt,temperture,counter);
    std::cout << "  [after]" << std::endl;
    for(int i = 0; i < num_particles; ++i){
      std::cout << "    i : " << i << " z[i] : " << z[i] << std::endl;
      std::cout << "    target E : " 
                << target_potential(i,z) 
                << std::endl;
    }
  }
  ASSERT_FALSE(z == z_previous);
} 

TEST_F(MetropolisEnsemblerTest, ThermodynamicCheck){
  std::size_t seed = 42;
  std::mt19937 mt(seed);
  int num_iteration = 1E+4;
  int dummy;
  int cooling_time = 100;
  double mean = 0;
  double var  = 0.0;
  Potential target_potential_temp = target_potential;

  Potential total_energy_func = 
    [target_potential_temp](int dummy, std::vector<double> const & z){
      int num_particles = z.size();
      double ene_t = 0; 
      for(int i = 0 ; i < num_particles; ++i) 
        ene_t += target_potential_temp(i,z);
      return ene_t/2;
    };

  std::vector<double> z_temp(num_particles,1);
  std::cout << "test total  E : " << total_energy_func(dummy,z_temp) 
            << " ,test target E : " << target_potential_temp(0,z_temp) 
                                    << std::endl;

  ensemble::Metropolis<Potential, Proposer> 
                    metropolis(target_potential,proposer);
  for(int cooling = 0; cooling < 10 * cooling_time;)
    metropolis.montecalro(z,mt,temperture,cooling);

  double ene_tt = 0;
  for(int counter = 0; counter < num_iteration; ++counter){
    for(int cooling = 0; cooling < cooling_time;)
      metropolis.montecalro(z,mt,temperture,cooling);
    ene_tt = total_energy_func(dummy, z);
    mean += ene_tt;
    var += ene_tt * ene_tt;
  }
  mean /= num_iteration;
  var /= num_iteration-1;
  var -= mean * mean * num_iteration / (num_iteration-1);
  var /= num_iteration;
  std::cout << "Total energy : " << mean 
            << " +/- " << 3 * std::pow(var,0.5) << std::endl;
}

}// end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
