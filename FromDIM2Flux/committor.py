"""Compute the committor function on milestones.

This program computes the committor function on milestones. The
computation is done by solving the committor equation for Markov
chains using sparse linear algebra.

The program uses two types of input files:

1. A matrix file containing a transition matrix sampled from an
equilibrium simulation (i.e., one where no boundary conditions have
been applied). This file is stored using the same file format used by
stationary-vectors and stationary-vectors.py.

2. A file containing indices (starting from 1) corresponding to a
special set of milestones. Each index should appear on a line by
itself. This is used to convey the identity of the reactant and
product milestones to the program.

For example, a reactant file containing the following two lines:

  1
  2

specifies that there are two milestones belonging to the reactant
state. Namely, the first two milestones (1 and 2).

"""

import sys

import numpy as np
import scipy.sparse
import scipy.sparse.linalg


__author__ = 'Juan M. Bello-Rivas'
__email__ = 'jmb@ices.utexas.edu'


def committor(K: scipy.sparse.csr_matrix,
              reactants: np.array,
              products: np.array) -> np.array:
    """Compute the committor function on the milestones.

    Parameters
    ----------
    K : scipy.sparse.csr_matrix
        Transition matrix coming from an equilibrium simulation (i.e.,
        a simulation in which no boundary conditions have been
        imposed).

    reactants : np.array
        Vector of boolean values identifying the indices corresponding
        to the milestones that belong to the reactant state.

    products : np.array
        Vector of boolean values identifying the indices corresponding
        to the milestones that belong to the product state.

    Returns
    -------
    committor_vector : np.array
        Vector of values of the committor function at the milestone
        corresponding to each index.

    """
    assert K.shape[0] == K.shape[1]
    M = K.shape[0]

    # Construct auxiliary matrix and solve for the values of the
    # committor function.
    reactant_indices = set(np.flatnonzero(reactants))
    product_indices = set(np.flatnonzero(products))
    reactant_and_product_indices = reactant_indices.union(product_indices)

    valid_indices = set(range(M)) - reactant_and_product_indices
    N = len(valid_indices)

    index_mapping = {i: k for k, i in enumerate(valid_indices)}

    A = scipy.sparse.lil_matrix((N, N))
    b = np.zeros(N)

    for i in valid_indices:
        ii = index_mapping[i]
        for j in range(M):
            Iij = 1 if i == j else 0
            Kij = K[i, j]

            if j in valid_indices:
                jj = index_mapping[j]
                A[ii, jj] = Iij - Kij
            elif j in product_indices:
                b[ii] = Kij - Iij

    c = scipy.sparse.linalg.spsolve(A.tocsr(), b)

    # Convert the solution into a vector whose entries are the values
    # of the committor function for the corresponding milestone.
    valid = np.fromiter((i in valid_indices for i in range(M)),
                        dtype=np.bool)

    committor_vector = np.zeros(M)
    committor_vector[valid] = c
    committor_vector[products] = 1.0

    return committor_vector


def load_matrix(filename: str) -> scipy.sparse.csr_matrix:
    """Load transition matrix from file.

    """
    def convert(x):
        return int(x[0])-1, int(x[1])-1, float(x[2])

    with open(filename) as f:
        M = int(f.readline().strip())
        lines = (line.split() for line in f.readlines())

    items = [np.array(convert(x)) for x in lines]
    ijv = np.array(items)

    K = scipy.sparse.csr_matrix((ijv[:, 2], (ijv[:, 0], ijv[:, 1])),
                                shape=(M, M))

    return K


def load_reactants_or_products(filename, M):
    """Load vector of reactant/product indices from file.

    This vector file contains indices corresponding to either a
    reactant or a product milestone, each on one line.

    """
    indices = np.loadtxt(filename, dtype=np.int)
    if not indices.shape:
        indices = {int(indices-1)}
    else:
        indices = set(indices-1)
    return np.fromiter((i in indices for i in range(M)), dtype=np.bool)


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Usage:', sys.argv[0], 'transition-matrix '
              'reactant-file product-file output-file',
              file=sys.stderr)
        sys.exit(-1)

    K = load_matrix(sys.argv[1])
    M = K.shape[0]
    reactants = load_reactants_or_products(sys.argv[2], M)
    products = load_reactants_or_products(sys.argv[3], M)

    committor_vector = committor(K, reactants, products)
    np.savetxt(sys.argv[4], committor_vector)
