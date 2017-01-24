#pragma once

#include <sys/types.h>

#include <map>
#include <vector>
#include <string>

namespace linear_algebra {

/**
 * Dense vector class.
 */
class vector {
 public:
  vector(unsigned dim);

  void zeros();
  void ones();

  double& operator[](unsigned k);
  double operator[](unsigned k) const;

  int save(std::string const& filename) const;

  inline size_t size() const { return dim; }

  double const* memptr() const { return &values[0]; }
  double* memptr() { return &values[0]; }

 private:
  const unsigned dim;
  std::vector<double> values;
};

/**
 * Sparse square matrix class.
 */
class matrix {
 public:
  matrix();

  int load(std::string const& filename);

  inline size_t size() const { return dim; }

  vector operator*(vector const& vec) const;
  void matrix_vector_multiply(double const* input, double* output) const;

 private:
  unsigned dim;
  std::map<int, std::vector<std::pair<int, double>>> entries;
};

int eigenpair(matrix const& A, vector& x, double* ritz_value = nullptr);

}
