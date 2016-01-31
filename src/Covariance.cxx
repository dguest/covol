#include "covol/Covariance.hh"
#include "H5Cpp.h"

#include <cassert>

Covariance::Covariance(const std::vector<CovVar>& vars):
  Covariance(vars.size())
{
  m_vars = vars;
}
Covariance::Covariance(const std::vector<std::string>& vars):
  Covariance(vars.size())
{
  for (const auto& name: vars) {
    m_vars.push_back( {name, ""} );
  }
}
// private constructor
Covariance::Covariance(size_t size):
  m_mean(Eigen::VectorXd::Zero(size)),
  m_comoment(Eigen::MatrixXd::Zero(size, size)),
  m_entries(0)
{
  assert(size > 0);
}

void Covariance::fill(const std::map<std::string, double>& vars, double wt) {
  if (wt == 0) return;
  using namespace Eigen;
  VectorXd invec(m_mean.size());
  for (int iii = 0; iii < m_mean.size(); iii++) {
    invec(iii) = vars.at(m_vars.at(iii).name);
  }
  double i = m_entries;
  VectorXd del = (invec - m_mean) / (i + 1);
  MatrixXd mat_del = i * del * del.transpose() - m_comoment / (i + 1);

  // update
  m_mean += del;
  m_comoment += mat_del;
  m_entries += wt;
}

Eigen::MatrixXd Covariance::getMatrix() const {
  // TODO: should I provide the unbiased covariance n / (n - 1)?
  return m_comoment;
}

std::ostream& operator<<(std::ostream& out, const Covariance& var) {
  out << "# Variables:\n";
  for (const auto& var: var.m_vars) {
    const auto& name = var.name;
    if (name.find(' ') != std::string::npos) {
      out << "'" <<  name << "' ";
    } else {
      out << name << " ";
    }
  }
  out << "\n";
  out << "# Cov Matrix:\n";
  out << var.getMatrix();
  return out;
}

namespace {
  struct H5Variable {
    const char* name;
    double mean;
    const char* units;
  };
  H5::CompType get_var_type() {
    // define types
    auto stype = H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);
    stype.setCset(H5T_CSET_UTF8);
    auto dtype = H5::PredType::NATIVE_DOUBLE;

    // define compound type
    H5::CompType type(sizeof(H5Variable));
    type.insertMember("name", offsetof(H5Variable, name), stype);
    type.insertMember("mean", offsetof(H5Variable, mean), dtype);
    type.insertMember("units", offsetof(H5Variable, units), stype);
    return type;
  }
  void add_variable_attributes(H5::DataSet& targ,
                               const std::vector<CovVar>& vars,
                               const Eigen::VectorXd& means) {
    const size_t size = vars.size();
    std::vector<H5Variable> h5_vars;
    for (int iii = 0; iii < size; iii++) {
      const auto& var = vars.at(iii);
      H5Variable hvar {var.name.c_str(), means(iii), var.units.c_str() };
      h5_vars.push_back(hvar);
    }
    std::vector<hsize_t> dim {size};
    H5::DataSpace space(1, dim.data());
    H5::CompType type = get_var_type();
    targ.createAttribute("variables", type, space).write(
      type, h5_vars.data());
  }
  void add_double_attribute(H5::DataSet& targ,
                            double val,
                            const std::string& name) {
    const auto type = H5::PredType::NATIVE_DOUBLE;
    targ.createAttribute(name, type, H5S_SCALAR).write(type, &val);
  }
}

void Covariance::write_to(H5::CommonFG& file,
                          const std::string& name,
                          int deflate) const {
  // setup dataspace
  size_t size = m_vars.size();
  std::vector<hsize_t> dims(2, size);
  H5::DataSpace ds(2, dims.data());

  // add properties
  H5::DSetCreatPropList params;
  params.setChunk(2, dims.data());
  params.setDeflate(deflate);

  // write the file
  const auto h5type = H5::PredType::NATIVE_DOUBLE;
  auto dataset = file.createDataSet(name, h5type, ds, params);
  dataset.write(getMatrix().data(), h5type);
  add_variable_attributes(dataset, m_vars, m_mean);
  add_double_attribute(dataset, m_entries, "sumwt");
}
