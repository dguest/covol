#ifndef COVARIANCE_HH
#define COVARIANCE_HH

#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

namespace H5 {
  class CommonFG;
}

class Covariance
{
public:
  Covariance(const std::vector<std::string>& variables);
  void fill(const std::map<std::string, double>& variables,
            double weight = 1.0);
  void write_to(H5::CommonFG& file,
                const std::string& name,
                int deflate = 7) const;
private:
  Eigen::MatrixXd getMatrix() const;
  std::vector<std::string> m_var_names;
  std::vector<std::string> m_var_units;
  double m_entries;
  Eigen::VectorXd m_mean;
  Eigen::MatrixXd m_comoment;
  friend std::ostream& operator<<(std::ostream&, const Covariance&);
};

std::ostream& operator<<(std::ostream&, const Covariance&);

#endif
