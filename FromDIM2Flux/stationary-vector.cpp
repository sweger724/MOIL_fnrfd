#include <cstdlib>

#include <string>
#include <vector>
#include <iostream>

#include "linear-algebra.h"
#include "timer.h"

using namespace std;

const std::string citation =
"****************************************************************************\n"
"This program uses a formula for computing the mean first passage time (MFPT)\n"
"that was introduced in:\n\n"
"  Bello-Rivas, J. M., & Elber, R. (2015). Exact milestoning. The Journal\n"
"  of Chemical Physics, 142(9), 094102. doi:10.1063/1.4913399\n\n"
"and further described in:\n\n"
"  Bello-Rivas, J. M., & Elber, R. (2015). Simulations of thermodynamics\n"
"  and kinetics on rough energy landscapes with milestoning. Journal of\n"
"  Computational Chemistry, doi:10.1002/jcc.24039.\n\n"
"If you use this program for computing the MFPT in a publication, please cite\n"
"the papers above.\n"
"****************************************************************************\n";

static double q_at_product(const std::vector<unsigned int>& product_indices,
                           const linear_algebra::vector& q);

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cerr << "Usage: " << argv[0] << " matrix-file lag-times-file stationary-vector-file\n";
    return EXIT_FAILURE;
  }

  timer clock;

  int status;

  linear_algebra::matrix m;

  clog << citation << "Reading matrix file \"" << argv[1] << "\"..." << endl;

  clock.tic();
  status = m.load(argv[1]);
  clock.toc();

  if (status != 0) {
    cerr << "Unable to load matrix into memory.\n";
    return EXIT_FAILURE;
  }

  clog << "Done reading matrix after " << clock.seconds() << " seconds.\n"
       << endl;

  linear_algebra::vector t;
  std::vector<unsigned> product_indices;
  status = t.load(argv[2], product_indices);
  if (status != 0 || product_indices.size() == 0) {
    cerr << "Something went wrong while loading the vector of lag times.\n";
    return EXIT_FAILURE;
  }

  linear_algebra::vector q(m.size());
  q.ones();

  if (linear_algebra::eigenpair(m, q) == 0)
    q.save(argv[3]);

  const double q_product = q_at_product(product_indices, q);
  const double mean_first_passage_time = q * t / q_product;

  cout << "Mean first passage time: " << mean_first_passage_time << endl;

  return EXIT_SUCCESS;
}


double q_at_product(const std::vector<unsigned int>& product_indices,
                    const linear_algebra::vector& q) {
  double sum = 0.0;

  for (auto i : product_indices)
    sum += q[i];

  return sum;
}
