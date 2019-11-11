#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>

#include"WangLandau.hpp"

double E_calc(const std::vector<double> z){
  double E = 0.0;
  double s_total_x = 0.0;
  double s_total_y = 0.0;
  for(int i = 0 ; i < z.size() ; ++i){ 
    s_total_x += std::cos(z[i]);
    s_total_y += std::sin(z[i]);
  }
  E = std::sqrt((s_total_x * s_total_x + s_total_y * s_total_y))/(double)z.size();
  std::cout << " E : " << E << std::endl;
  return E;
}

void Proposer(std::vector<double>& z, std::mt19937 & mt){
  std::uniform_real_distribution<> uniform_rand(0.0,2.0*M_PI);
  for(int i = 0 ; i < z.size(); ++i) z[i] = uniform_rand(mt);
}

int main() {

  std::cout << "[test " << metadynamics::WangLandau::name() << "]\n";
  std::size_t seed = 1234;
  std::mt19937 mt(seed);

  //declare class 
  unsigned int dimension = 100;
  double bin = 0.01;
  double E_max = 1.0;
  double E_min = 0.0;
  unsigned int N_bin = int((E_max - E_min) / bin);
  int max_iteration = 100;
  std::cout << "[declare class]" << std::endl;
  metadynamics::WangLandau wang_landau(dimension, N_bin, E_max, E_min);
  //end declare class

  std::vector<double> density(N_bin),entropy(N_bin),bins(N_bin);
  std::cout << "[start sampling]" << std::endl;
  wang_landau.montecalro(E_calc,Proposer,max_iteration,mt);
  std::cout << "[end sampling]" << std::endl;
  bins = wang_landau.bins();
  density = wang_landau.density();
  entropy = wang_landau.entropy();

  std::cout << "[result put]" << std::endl;
  std::ofstream ofs("result_density.dat",std::ios::out);
  std::ofstream ofsr("result_entropy.dat",std::ios::out);
  ofs << "#output bins density exact" << std::endl;
  ofsr << "#output bins entropy exact" << std::endl;
  for(int i = 0 ; i < N_bin ; ++i) ;
  auto exact = [](double x) {return 1;};
  double c_norm = 0.0;
  for(int i = 0 ; i < N_bin ; ++i){
    ofs << bins[i] << " " << density[i] << " " << exact(bins[i]) << std::endl;
    c_norm += density[i] * (bins[0] - bins[1]);
    ofsr << bins[i] << " " << entropy[i] << " " << std::log(exact(bins[i])) << std::endl;
  }
  std::cout << " c_norm : " << c_norm << std::endl;
  std::cout << "[finish test]" << std::endl;

  return 0;
}
