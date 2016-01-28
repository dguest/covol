#include "covol/Covariance.hh"
#include <cstdlib>
#include <iostream>
#include <random>

#include "H5Cpp.h"

int main(int argc, char* argv[]) {
  int num = 10;
  if (argc > 1) num = atoi(argv[1]);

  std::default_random_engine random_gen;
  std::normal_distribution<double> norm(0,1);

  std::vector<std::string> vars {
    "long_long_long long",
      "pt", "eta", "1", "2", "3", "4", "5"};
  Covariance cov(vars);
  for (int iii = 0; iii < num; iii++) {
    double rand = norm(random_gen) * 1;
    cov.fill( {
        {"long_long_long long", norm(random_gen)}, 
        {"pt", norm(random_gen) + rand},
        {"1", norm(random_gen) + rand},
        {"2", norm(random_gen) + rand},
        {"3", norm(random_gen) + rand},
        {"4", norm(random_gen) + rand},
        {"5", norm(random_gen) + rand},
        {"eta", norm(random_gen) - rand} } );
  }
  std::cout << cov << std::endl;

  H5::H5File out_file("test.h5", H5F_ACC_TRUNC);
  cov.write_to(out_file, "cov_matrix");
  return 0;
}
