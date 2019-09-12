#include <gtest/gtest.h>
#include <lattice/chain.hpp>

TEST(ChainTest, BasicTest){
  // set lattice 
  int system_size = 10;
  lattice::Chain lattice(system_size);
  int num_particles = lattice.set_num_particles(system_size); 
  int num_adjacent = lattice.number_adjacent();
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(num_adjacent));;
  lattice.create_table(pair_table);
  EXPECT_EQ(num_particles, system_size);
  EXPECT_EQ(num_adjacent,2);
  EXPECT_EQ(pair_table[0][0],num_particles-1);
  EXPECT_EQ(pair_table[0][1],1);
  EXPECT_EQ(pair_table[num_particles-1][0],num_particles-2);
  EXPECT_EQ(pair_table[num_particles-1][1],0);
}
