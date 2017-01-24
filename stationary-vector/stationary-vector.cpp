#include <cstdlib>

#include <iostream>

#include "linear-algebra.h"
#include "timer.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " matrix-file stationary-vector-file\n";
    return EXIT_FAILURE;
  }

  timer clock;

  int status;

  linear_algebra::matrix m;

  clog << "Reading matrix file \"" << argv[1] << "\"..." << endl;

  clock.tic();
  status = m.load(argv[1]);
  clock.toc();

  if (status != 0) {
    cerr << "Unable to load matrix into memory.\n";
    return EXIT_FAILURE;
  }

  clog << "Done reading matrix after " << clock.seconds() << " seconds.\n"
       << endl;

  linear_algebra::vector v(m.size());
  v.ones();

  if (linear_algebra::eigenpair(m, v) == 0)
    v.save(argv[2]);

  return EXIT_SUCCESS;
}
