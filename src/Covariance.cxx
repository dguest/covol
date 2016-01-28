#include "Covariance.hh"

// #include <iostream>

Covariance::Covariance(const std::vector<std::string>& vars):
  m_var_names(vars),
  m_entries(0),
  m_mean(Eigen::VectorXd::Zero(vars.size())),
  m_comoment(Eigen::MatrixXd::Zero(vars.size(), vars.size()))
{
}

void Covariance::fill(const std::map<std::string, double>& vars) {
  using namespace Eigen;
  VectorXd invec(m_mean.size());
  for (int iii = 0; iii < m_mean.size(); iii++) {
    invec(iii) = vars.at(m_var_names.at(iii));
  }
  double i = m_entries;
  VectorXd del = (invec - m_mean) / (i + 1);
  MatrixXd mat_del = i * del * del.transpose() - m_comoment / (i + 1);

  // update
  m_mean += del;
  m_comoment += mat_del;
  m_entries++;
}

Eigen::MatrixXd Covariance::getMatrix() const {
  // TODO: should I provide the unbiased covariance n / (n - 1)?
  return m_comoment;
}

std::vector<std::string> Covariance::getVariables() const {
  return m_var_names;
}
