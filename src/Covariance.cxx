#include "covol/Covariance.hh"
#include "H5Cpp.h"

Covariance::Covariance(const std::vector<std::string>& vars):
  m_var_names(vars),
  m_entries(0),
  m_mean(Eigen::VectorXd::Zero(vars.size())),
  m_comoment(Eigen::MatrixXd::Zero(vars.size(), vars.size()))
{
}

void Covariance::fill(const std::map<std::string, double>& vars, double wt) {
  using namespace Eigen;
  VectorXd invec(m_mean.size());
  for (int iii = 0; iii < m_mean.size(); iii++) {
    invec(iii) = vars.at(m_var_names.at(iii));
  }
  double i = m_entries;
  VectorXd x = invec;
  MatrixXd mat = m_comoment;
  VectorXd del = (x - m_mean) * wt / (i + wt);
  MatrixXd mat_del = i / wt * del * del.transpose() - mat * wt / (i + wt);

  // update
  m_mean += del;
  m_comoment += mat_del;
  m_entries += wt;
}

Eigen::MatrixXd Covariance::getMatrix() const {
  // TODO: should I provide the unbiased covariance n / (n - 1)?
  // FIXME: do I need this?
  return m_comoment;
}

std::ostream& operator<<(std::ostream& out, const Covariance& var) {
  out << "# Variables:\n";
  for (const auto& name: var.m_var_names) {
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
    // TODO: add units?
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
    return type;
  }
  void add_variable_attributes(H5::DataSet& targ,
                               const std::vector<std::string>& names,
                               const Eigen::VectorXd& means) {
    const size_t size = names.size();
    std::vector<H5Variable> h5_vars;
    for (int iii = 0; iii < size; iii++) {
      H5Variable var {names.at(iii).c_str(), means(iii) };
      h5_vars.push_back(var);
    }
    std::vector<hsize_t> dim {size};
    H5::DataSpace space(1, dim.data());
    H5::CompType type = get_var_type();
    targ.createAttribute("variables", type, space).write(
      type, h5_vars.data());
  }
}

void Covariance::write_to(H5::CommonFG& file,
                          const std::string& name,
                          int deflate) const {
  // setup dataspace
  size_t size = m_var_names.size();
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
  add_variable_attributes(dataset, m_var_names, m_mean);
}
