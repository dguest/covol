#ifndef COVARIANCE_HH
#define COVARIANCE_HH

#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

namespace H5 {
  class CommonFG;
}

struct CovVar {
  std::string name;
  std::string units;
};

class Covariance
{
public:
  Covariance(const std::vector<CovVar>& variables);
  // convenience constructor (no units)
  Covariance(const std::vector<std::string>& variables);
  void fill(const std::map<std::string, double>& variables,
            double weight = 1.0);
  void write_to(H5::CommonFG& file,
                const std::string& name,
                int deflate = 7) const;
private:
  Covariance(size_t);
  Eigen::MatrixXd getMatrix() const;
  Eigen::VectorXd m_mean;
  Eigen::MatrixXd m_comoment;
  double m_entries;
  std::vector<CovVar> m_vars;
  friend std::ostream& operator<<(std::ostream&, const Covariance&);
};

std::ostream& operator<<(std::ostream&, const Covariance&);

#endif
