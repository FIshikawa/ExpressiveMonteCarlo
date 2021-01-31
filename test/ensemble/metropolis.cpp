#include <gtest/gtest.h>
#include <limits>
#include <lattice/chain.hpp>
#include <ensemble/metropolis.hpp>


using Potential = std::function<double(int,std::vector<double> const &)>;
using Proposer = std::function<std::pair<int,double>(std::vector<double> const &,std::mt19937 &)>;

namespace {
class MetropolisEnsemblerTest : public ::testing::Test {
protected:
  virtual void SetUp(){
    // set lattice 
    int system_size = 10;
    lattice::Chain lattice(system_size);
    num_particles = lattice.set_num_particles(system_size); 
    num_adjacent = lattice.number_adjacent();
    pair_table = std::vector<std::vector<int> >
                    (num_particles,std::vector<int>(num_adjacent));
    lattice.create_table(pair_table);

    // set constants
    J = 1.0;
    precision = std::numeric_limits<double>::epsilon();
    temperture = 1.0;

    // set vector
    z.resize(num_particles);
    for(int i = 0; i < num_particles; ++i) z[i] = -1;

    // set functions
    // set target potential
    int num_particles_temp = num_particles;
    proposer = 
      [num_particles_temp](std::vector<double> const & z, std::mt19937 & mt){
        std::uniform_real_distribution<> index_dist(0,num_particles_temp-1);
        int target = index_dist(mt);
        double new_value = ((z[target] < 0) ? 1 : -1);
        std::pair<int,double> updates = std::make_pair(target, new_value);
        return updates;
      };

    // set proposer 
    int num_adjacent_temp = num_adjacent;
    double J_temp = J;
    std::vector<std::vector<int> > pair_table_temp = pair_table;
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
  double J,precision,temperture;
  std::vector<std::vector<int> > pair_table;
};

TEST_F(MetropolisEnsemblerTest, BasicTest){
  int counter = 0;
  std::size_t seed = 1234;
  std::mt19937 mt(seed);

  ensemble::Metropolis<Potential, Proposer> metropolis(target_potential,proposer);
  int num_iteration = 2;
  while(counter < num_iteration){
    std::cout << "<< counter : " << counter << " >>" << std::endl;
    std::cout << "  [before]" << std::endl;
    for(int i = 0; i < num_particles; ++i)
      std::cout << "    i : " << i << " z[i] : " << z[i] << std::endl;
    metropolis.montecalro(z,mt,temperture,counter);
    std::cout << "  [after]" << std::endl;
    for(int i = 0; i < num_particles; ++i)
      std::cout << "    i : " << i << " z[i] : " << z[i] << std::endl;
  }
  ASSERT_EQ(counter,num_iteration);
} 

TEST_F(MetropolisEnsemblerTest, ThermodynamicCheck){
  std::size_t seed = 1234;
  std::mt19937 mt(seed);
  int num_iteration = 1E+4;
  int dummy;
  int cooling_time = 10;
  double mean = 0.0;
  double var  = 0.0;
  Potential target_potential_temp = target_potential;

  Potential total_energy_func = 
    [target_potential_temp](int dummy, std::vector<double> const & z){
      int num_particles = z.size();
      double ene = 0; 
      for(int i = 0 ; i < num_particles; ++i) 
        ene += target_potential_temp(i,z);
      ene /= 2;
      return ene;
    };

  ensemble::Metropolis<Potential, Proposer> metropolis(target_potential,proposer);
  int previous;
  for(int counter = 0; counter < num_iteration;){
    for(int cooling = 0; cooling < cooling_time;)
      metropolis.montecalro(z,mt,temperture,cooling);
    previous = counter;
    metropolis.montecalro(z,mt,temperture,counter);
    if(previous < counter){
      double ene = total_energy_func(dummy, z);
      mean += ene;
      var += ene * ene;
    }
  }
  mean /= num_iteration;
  var /= num_iteration;
  std::cout << "Total energy : " << total_energy << std::endl;
}

}// end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
