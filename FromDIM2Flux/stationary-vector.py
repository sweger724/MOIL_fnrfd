#!/usr/bin/env python


import sys
import math
import logging
import warnings
from typing import Tuple

import numpy as np
import scipy.sparse
import scipy.sparse.linalg as linalg
from scipy.sparse.linalg.eigen.arpack.arpack import ArpackError


EPSILON = sys.float_info.epsilon

SHIFT = 1.0 + EPSILON

DOMINANT_EIGENVALUE_TOLERANCE = 1e-3

SPECTRAL_GAP_TOLERANCE = 1e-9


def dominant_eigenvector(transition_matrix: scipy.sparse.csr_matrix,
                         normalize: bool = True) \
        -> Tuple[float, np.array]:
    '''Compute dominant eigenvector of a row-stochastic matrix.

    Parameters
    ----------

    transition_matrix
        A row-stochastic matrix in Compressed Sparse Row format.

    normalize
        Whether to normalize the resulting eigenvector so that its
        1-norm is equal to one.

    Returns
    -------

    eigenvalue
        Dominant eigenvalue.

    eigenvector
        The eigenvector corresponding to the dominant eigenvalue.

    Raises
    ------

    ValueError
        If something goes wrong during the computation.

    '''
    # assert sparse.isspmatrix_csr(transition_matrix)

    try:
        v = np.ones(transition_matrix.shape[0])
        eigenvalues, eigenvectors = linalg.eigs(transition_matrix,
                                                k=2, sigma=SHIFT,
                                                which='LM', v0=v)
    except ArpackError as exc:
        logging.debug('The sparse eigensolver failed: "{}" This '
                      'sometimes occurs for matrices with less '
                      'than 10 non-empty rows or columns.'.format(exc))
        eigenvalues, eigenvectors = dense_eigensolver(transition_matrix)

    if eigenvalues[0].real > eigenvalues[1].real:
        largest, next_largest = 0, 1
    else:
        largest, next_largest = 1, 0

    if not math.isclose(eigenvalues[largest].imag, 0.0,
                        abs_tol=1e1 * EPSILON):
        warnings.warn('Dominant eigenvalue {} is not a real number'
                      .format(eigenvalues[-1]), RuntimeWarning)

    spectral_gap = np.abs(eigenvalues[largest] -
                          eigenvalues[next_largest])
    if math.isclose(spectral_gap, 0.0, abs_tol=SPECTRAL_GAP_TOLERANCE):
        warnings.warn('Spectral gap {} is smaller than the prescribed '
                      'tolerance of {}'.format(spectral_gap,
                                               SPECTRAL_GAP_TOLERANCE),
                      RuntimeWarning)

    eigenvalue = eigenvalues[largest].real
    eigenvector = eigenvectors[:, largest].real

    if not math.isclose(eigenvalue, 1.0,
                        abs_tol=DOMINANT_EIGENVALUE_TOLERANCE):
        warnings.warn('Dominant eigenvalue {} is not close enough to 1'
                      .format(eigenvalue), RuntimeWarning)

    if normalize is True:
        sgn = np.sign(eigenvector.sum())
        eigenvector = sgn * eigenvector / np.linalg.norm(eigenvector, 1)

    return eigenvalue, ensure_nonnegative(eigenvector)


def ensure_nonnegative(vector: np.array,
                       tolerance: float = EPSILON) \
        -> np.array:
    '''Ensure all entries of vector are nonnegative.

    If an entry is negative but smaller in magnitude than the value of
    the tolerance parameter, then we replace it by its absolute value.

    '''
    assert -tolerance < np.min(vector), \
        ('Element with value {} surpasses'
         'tolerance {}'.format(np.min(vector), -tolerance))
    return np.abs(vector)


def rayleigh_quotient(matrix: scipy.sparse.csr_matrix,
                      vector: np.array) -> float:
    '''Compute the Rayleigh quotient of a matrix and a vector.

    '''
    return matrix.dot(vector).T.dot(vector) / vector.dot(vector)


def dense_eigensolver(matrix: scipy.sparse.csr_matrix) \
        -> Tuple[np.array, np.array]:
    '''Compute the two largest eigenvalues and the corresponding
    eigenvectors using dense linear algebra.

    Parameters
    ----------
    matrix
        Sparse MxM matrix in compressed sparse row or column format.

    Returns
    -------
    eigenvalues
        The first two largest eigenvalues
    eigenvectors
        Eigenvectors corresponding to the two largest eigenvalues.
        Use the last index to access each eigenvector.

    '''
    # assert scipy.sparse.isspmatrix_csr(matrix), 'Invalid sparse matrix'

    # Convert to a dense matrix.
    K = matrix
    row, col = K.nonzero()
    indices = set(row).union(set(col))
    M = len(indices)

    map_forward = {}
    map_backward = {}
    for i, idx in enumerate(indices):
        map_forward[idx] = i
        map_backward[i] = idx

    KK = np.zeros((M, M), dtype=np.float64)
    for i, j in zip(row, col):
        ii, jj = map_forward[i], map_forward[j]
        KK[ii, jj] = K[i, j]

    # Solve eigenproblem using dense linear algebra.
    ew, ev = np.linalg.eig(KK)
    permutation = np.argsort(ew)
    ew = ew[permutation]
    ev = ev[:, permutation]

    eigenvalues = ew[-2:]

    # Convert the dominant eigenvector back to the original form.
    eigenvectors = np.zeros((K.shape[0], 2), dtype=K.dtype)
    for j in [-2, -1]:
        for ii in range(M):
            i = map_backward[ii]
            eigenvectors[i, j] = ev[ii, j].real

    return eigenvalues, eigenvectors


CITATION = '''\
*****************************************************************************
This program uses a formula for computing the mean first passage time (MFPT)
that was introduced in:

  Bello-Rivas, J. M., & Elber, R. (2015). Exact milestoning. The Journal
  of Chemical Physics, 142(9), 094102. doi:10.1063/1.4913399

and further described in:

  Bello-Rivas, J. M., & Elber, R. (2015). Simulations of thermodynamics
  and kinetics on rough energy landscapes with milestoning. Journal of
  Computational Chemistry, doi:10.1002/jcc.24039.

If you use this program to produce results for a publication, please cite the
papers above.
*****************************************************************************

'''


def load_matrix(filename):
    '''Load matrix from file.'''
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


def load_vector(filename, M):
    '''Load vector of lag times from file.'''
    t = np.zeros(M)

    tt = np.loadtxt(filename)
    indices = tt[:, 0].astype(int)-1
    values = tt[:, 1]

    t[indices] = values
    product_indices = (t == 0).nonzero()[0]

    return t, product_indices


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: {} transition-matrix local-lag-times '
              'dominant-eigenvector'.format(sys.argv[0]),
              file=sys.stderr)
        sys.exit(-1)

    # Load transition matrix.
    try:
        transition_matrix = load_matrix(sys.argv[1])
    except IOError as exc:
        print('Error reading transition matrix file:', exc)
        sys.exit(-1)

    l, q = dominant_eigenvector(transition_matrix.T)

    # Load vector of lag times.
    M = transition_matrix.shape[0]
    t, product_indices = load_vector(sys.argv[2], M)

    # Compute MFPT.
    sum_products = q[product_indices].sum()
    mfpt = q.dot(t) / sum_products
    print('{}The estimated MFPT is equal to {}'.format(CITATION, mfpt))

    # Save stationary vector.
    with open(sys.argv[3], 'w') as f:
        for i in range(q.shape[0]):
            print(q[i], file=f)
