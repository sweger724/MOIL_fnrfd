#include <cassert>

#include <fstream>
#include <algorithm>

#include "linear-algebra.h"

namespace linear_algebra {

vector::vector()
    : dim(0),
      values(dim, 0.0) {
}

vector::vector(unsigned dim_)
    : dim(dim_),
      values(dim, 0.0) {
}

void vector::zeros() {
  std::fill(values.begin(), values.end(), 0.0);
}

void vector::ones() {
  std::fill(values.begin(), values.end(), 1.0);
}

double& vector::operator[](unsigned k) {
  assert(k < dim);
  return values[k];
}

double vector::operator[](unsigned k) const {
  assert(k < dim);
  return values[k];
}

int vector::save(std::string const& filename) const {
  std::ofstream ofs(filename);
  if (!ofs.is_open())
    return -1;

  for (auto x : values) {
    ofs << x << "\n";
    if (!ofs.good())
      return -1;
  }

  return 0;
}

int vector::load(std::string const& filename, std::vector<unsigned>& product_indices) {
  std::ifstream ifs(filename);
  if (!ifs.is_open())
    return -1;

  values = std::vector<double>();

  while (ifs.good()) {
    int i;
    double v;

    ifs >> i >> v;

    if (ifs.eof())
      break;

    if (v == 0) product_indices.push_back(i - 1);

    values.push_back(v);
  }

  dim = values.size();

  return 0;
}

double vector::operator*(const vector& other) {
  assert(other.dim == dim);

  double result = 0.0;

  for (unsigned i = 0; i < dim; ++i)
    result += (*this)[i] * other[i];

  return result;
}

matrix::matrix() : dim(0) {}

int matrix::load(std::string const& filename) {
  std::ifstream ifs(filename);
  if (!ifs.is_open())
    return -1;

  int N;
  ifs >> N;
  if (N <= 0)
    return -1;

  dim = N;

  while (ifs.good()) {
    int i, j;
    double aij;

    ifs >> i >> j >> aij;

    if (ifs.eof())
      break;

    if (aij == 0.0)
      continue;

    // Be careful: the matrix is transposed when stored in memory!
    entries[j-1].push_back(std::make_pair(i-1, aij));
  }

  return 0;
}

vector matrix::operator*(vector const& vec) const {
  assert(vec.size() == dim);

  vector w(dim);

  matrix_vector_multiply(vec.memptr(), w.memptr());

  return w;
}

void matrix::matrix_vector_multiply(double const* input, double* output) const {
  assert(input != nullptr);
  assert(output != nullptr);
  assert(dim > 0);

  std::fill_n(output, dim, 0.0);

  for (auto const& ij : entries) {
    const int i = ij.first;

    for (auto const& j_aij : ij.second) {
      const int j = j_aij.first;
      const double aij = j_aij.second;

      output[i] += aij * input[j];
    }
  }
}

}
