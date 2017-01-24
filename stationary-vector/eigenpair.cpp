#include <cmath>
#include <cstring>

#include <limits>
#include <vector>
#include <iostream>
#include <iomanip>

#include "arnoldi.h"
#include "linear-algebra.h"
#include "timer.h"

namespace linear_algebra {

const double epsilon = std::numeric_limits<double>::epsilon();

// Obtain dominant eigenpair.
int eigenpair(matrix const& A, vector& x0, double* ritz_value) {
  int n = A.size();
  int nev = 1;
  int ncv = 4;
  int maxn = n;
  int maxncv = ncv;
  int ldv = maxn;

  double tol = epsilon;

  std::vector<double> d(maxncv * 3);
  double* resid = x0.memptr();

  std::vector<double> v(ldv * maxncv);
  std::vector<double> workd(3 * maxn);
  std::vector<double> workev(3 * maxncv);
  std::vector<double> workl(3 * maxncv * maxncv + 6 * maxncv);

  std::vector<int> iparam(11);
  std::vector<int> ipntr(14);
  std::vector<char> select(maxncv);

  char bmat[] = "I";
  char which[] = "LM";

  int ido = 0;
  int lworkl= 3 * ncv * ncv + 6 *ncv;
  int info = 1;

  int ishfts = 1;
  int maxitr = 1e7;
  int mode1 = 1;

  iparam[0] = ishfts;
  iparam[2] = maxitr;
  iparam[6] = mode1;

  std::clog << "Running Arnoldi iteration with tolerance " << tol << "...\n";

  timer clock;
  clock.tic();

  for (int itr = 0; ido == 99 || itr < maxitr; ++itr) {
    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, &v[0], &ldv,
            &iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &info);

    if (info < 0) {
      std::cerr << "Error after call to dnaupd: info = " << info << std::endl;
      return -1;
    }

    if (ido != 1)
      break;

    A.matrix_vector_multiply(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
  }

  clock.toc();
  std::clog << "Done with Arnoldi iteration (" << clock.seconds() << " seconds)."
            << std::endl << std::endl;

  int ierr;
  int rvec = 1;
  char how_many = 'A';
  double sigmar = 0.0, sigmai = 0.0;

  std::clog << "Running postprocessing step...\n";
  clock.tic();

  dneupd_(&rvec, &how_many, &select[0], &d[0], &d[maxncv], &v[0], &ldv, &sigmar,
          &sigmai, &workev[0], bmat, &n, which, &nev, &tol, resid, &ncv, &v[0], &ldv,
          &iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &ierr);

  clock.toc();

  int nconv = iparam[4];

  std::clog << "Done with postprocessing step. Status: "
            << (ierr == 0 ? "OK" : "FAIL")
            << " (" << clock.seconds() << " seconds)." << "\n"
            << "Number of converged Ritz values: " << nconv << "\n"
            << "Number of implicit Arnoldi updates: " << iparam[2] << "\n"
            << "Number of matrix-vector products: " << iparam[8] << std::endl;

  if (ierr != 0) {
    std::clog << "ierr = " << ierr << std::endl;
    return -1;
  }

  // Normalize result.
  double norm1 = 0.0, sum = 0.0;
  for (int i = 0; i < maxn; ++i) {
    x0[i] = v[i + 0 * maxn];
    norm1 += std::fabs(x0[i]);
    sum += x0[i];
  }

  const double sign = std::signbit(sum) == 1 ? -1.0 : 1.0;
  for (int i = 0; i < maxn; ++i)
    x0[i] /= (norm1 * sign);
  
  // Save the largest Ritz value (= best approximation to the largest
  // eigenvalue).
  std::clog << "Largest Ritz value: " << d[0] << std::endl;
  if (ritz_value != nullptr)
    *ritz_value = d[0];

  return 0;
}

}
