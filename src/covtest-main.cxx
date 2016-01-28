#include "Covariance.hh"
#include <cstdlib>
#include <iostream>
#include <random>

#include "H5Cpp.h"

int main(int argc, char* argv[]) {
  int num = 10;
  if (argc > 1) num = atoi(argv[1]);

  std::default_random_engine random_gen;
  std::normal_distribution<double> norm(0,1);

  std::vector<std::string> vars {"pt", "eta"};
  Covariance cov(vars);
  for (int iii = 0; iii < num; iii++) {
    double rand = norm(random_gen) * 0.1;
    cov.fill( {
        {"pt", norm(random_gen) + rand},
        {"eta", norm(random_gen) - rand} } );
  }
  H5::H5File out_file("test.h5", H5F_ACC_TRUNC);
  cov.write_to(out_file, "cov_matrix");
  return 0;
}
