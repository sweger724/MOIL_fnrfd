/////////////////////////////////////////////////////////////////////////////
//
// edward.cpp -- Estimates the variance of the first passage times by
// resampling from the Markov chain determined by the transition
// matrix K and the vector of local MFPTs t.
//
/////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cassert>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>

#include "linear-algebra.h"

typedef std::vector<size_t> indices;

class edward {
 public:
  edward(linear_algebra::matrix& K_,
         linear_algebra::vector& t_,
         const indices& reactants, const indices& products)
      : K(K_), t(t_),
        reactant_indices(reactants), product_indices(products),
        p0(0, reactants.size() - 1), unif(0.0, 1.0) {
  }

  double sample(std::mt19937& rng);

 private:
  linear_algebra::matrix K;
  linear_algebra::vector t;
  const indices& reactant_indices, product_indices;
  std::uniform_int_distribution<size_t> p0;
  std::uniform_real_distribution<double> unif;
};

static void print_indices(std::ostream& ostr, indices idcs, std::string title) {
  ostr << title << " ";
  for (auto idx : idcs)
    ostr << idx+1 << " ";
  ostr << "\n";
}

int main(int argc, char* argv[]) {
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0] << " matrix-file local-mfpt-file"
              << " num-resamples random-seed"
              << " reactant-index-1 reactant-index-2...\n";
    return EXIT_FAILURE;
  }

  int arg = 0;
  const std::string matrix_filename(argv[++arg]);
  const std::string vector_filename(argv[++arg]);

  const size_t num_resamples = atof(argv[++arg]);
  if (num_resamples < 1) {
    std::cerr << "Error: invalid number of resamples.  "
              << "It should be at least 1.\n";
    return EXIT_FAILURE;
  }

  const unsigned long random_seed = strtoul(argv[++arg], nullptr, 10);

  std::clog << "Using random seed: " << random_seed << std::endl;

  std::cout << std::setprecision(9);
  std::clog << std::setprecision(9);

  indices reactant_indices;
  {
    for (int k = arg+1; k < argc; ++k)
      reactant_indices.push_back(atoi(argv[k])-1);
  }

  linear_algebra::matrix K;
  if (K.load(matrix_filename, false) != 0) {
    std::cerr << "Error: unable to load transition matrix `"
              << argv[1] << "'\n";
    return EXIT_FAILURE;
  }

  linear_algebra::vector t;
  indices product_indices;
  if (t.load(vector_filename, product_indices) != 0) {
    std::cerr << "Error: unable to load vector of local MFPTs `"
              << argv[2] << "'\n";
    return EXIT_FAILURE;
  }

  edward hopper(K, t, reactant_indices, product_indices);

  print_indices(std::clog, reactant_indices, "Reactant indices:");
  print_indices(std::clog, product_indices, "Product indices:");

  std::mt19937 rng(random_seed);

  std::vector<double> fpts;
  for (size_t n = 0; n < num_resamples; ++n)
    fpts.push_back(hopper.sample(rng));

  const double mfpt = std::accumulate(fpts.begin(), fpts.end(), 0.0) / double(fpts.size());
  double vfpt = 0.0;
  for (auto fpt : fpts)
    vfpt += (fpt - mfpt) * (fpt - mfpt);
  vfpt /= double(fpts.size() - 1);
  std::cout << "FPT mean: " << mfpt << " std.: " << std::sqrt(vfpt) << "\n";

  return EXIT_SUCCESS;
}

static inline bool reached_product(const indices& product_indices, size_t i) {
  return std::find(product_indices.begin(), product_indices.end(), i)
      != product_indices.end();
}

double edward::sample(std::mt19937& rng) {
  // Caveat: as it stands, the initial milestone is uniformly
  // distributed while the correct distribution should be based on the
  // relative (canonical) weights of the milestones belonging to the
  // reactant.
  const size_t start_idx = reactant_indices[p0(rng)];

  double fpt = 0.0;
  size_t i;
  ssize_t j;
  for (i = start_idx; ; i = j) {
    if (reached_product(product_indices, i))
      break;

    fpt += t[i];

    const double c = unif(rng);
    double sum = 0.0;
    j = -1;
    for (auto elem : K.get_row(i)) {
      j = std::get<0>(elem);
      auto val = std::get<1>(elem);
      sum += val;
      if (sum >= c)
        break;
    }

    assert(sum >= c && j >= 0);
  }

  return fpt;
}
