#include <gtest/gtest.h>
#include <limits>
#include <lattice/chain.hpp>
#include <ensemble/metropolis.hpp>

namespace {
class MetropolisEnsemblerTest : public ::testing::Test {
protected:
  virtual void SetUp(){
    // set lattice 
    int system_size = 10;
    lattice::Chain lattice(system_size);
    num_particles = lattice.set_num_particles(system_size); 
    ASSERT_EQ(num_particles, system_size);
    num_adjacent = lattice.number_adjacent();
    int num_adjacent_temp = num_adjacent;
    ASSERT_EQ(num_adjacent,2);
    pair_table = std::vector<std::vector<int> >(num_particles,std::vector<int>(num_adjacent));
    lattice.create_table(pair_table);
    ASSERT_EQ(pair_table[0][0],num_particles-1);
    ASSERT_EQ(pair_table[0][1],1);
    ASSERT_EQ(pair_table[num_particles-1][0],num_particles-2);
    ASSERT_EQ(pair_table[num_particles-1][1],0);

    // set constants
    J = 1.0;
    precision = std::numeric_limits<double>::epsilon();

    // set vector
    z.resize(num_particles);
    for(int i = 0; i < num_particles; ++i) z[i] = -1;
 
    // set proposer
    std::uniform_real_distribution<> index_dist(0,num_particles-1);
    proposer = [&index_dist](std::vector<double>z, std::mt19937 mt){
      int target = index_dist(mt);
      double new_value = (z[target] < 0) ? 1 : -1;
      return std::make_pair<int,double>(target, new_value);
    };

    // set potential energy
    double J_temp = J;
    std::vector<std::vector<int> > pair_table_temp = pair_table;
    target_potential = [&J_temp,&num_adjacent_temp,&pair_table_temp](int target, std::vector<double> z, double t){
      double ene = 0; 
      for(int i = 0 ; i < num_adjacent_temp; ++i) ene += J_temp*z[target]*z[pair_table_temp[target][i]];
      return ene;
    };

    // set random numbers
    std::size_t seed = 1234;
    std::mt19937 mt(seed);
  } 
  std::vector<double> z;
  int num_particles, num_adjacent;
  double J,precision;
  std::vector<std::vector<int> > pair_table;
  std::function<double(int,std::vector<double>,double)> target_potential;
  std::function<std::pair<int,double>(std::vector<double>, std::mt19937)> proposer; 
};

}// end namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
}
