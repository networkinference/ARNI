#!/usr/bin/env python
import numpy as np

def basis_expansion(X, K, TYPE, NODE):
    '''
     basis_expansion(X,K,TYPE,NODE) generates a multidimensional array of
     basis expansions evaluated on all points of a multivariate time series.

     Parameters
     ------------------
     X:    Matrix containing N time series of M time points.
     K:    Maximum order of the basis expansion.
     TYPE: Type of basis function employed. In this file, we only
           implement expansions up to pairwise interactions. Expansions
           availables are: polynomial (x_{j}^{k}), polynomial_diff
           ((x_{j}-x_{NODE})^{k}), power_series (x_{j}^{k1}*x_{i}^{k2}),
           fourier (sin(k*x_{j}) and cos(k*x_{j})), fourier_diff
           (sin(k*(x_{j}-x_{NODE})) and cos(k*(x_{j}-x_{NODE}))) and RBF (a model
           based on radial basis functions). These functions are shown in
           table I in the main manuscript.
     NODE: Unit on which we are performing the reconstruction. Zero indexed.

     Input type
     ------------------
     X:    double
     K:    integer
     TYPE: string
     NODE: integer

     Output
     ------------------
     Expansion: Multidimensional array of size [K+1,M,N] containing the
     evalation of all k=0,1,...,K basis functions for all M time points and
     all N possible incoming connections. For power_series, (K*K+2) basis
     functions are employed, and for fourier(_diff), 2*(K+1) are employed.

     Example
     ------------------
     basis_expansion(X,4,'power_series',5); generates a multidimensional array
     of size [18,M,N] containing the evaluation of the basis for all M time
     points and all N possible incoming connections.

     Accompanying material to "Model-free inference of direct interactions
     from nonlinear collective dynamics".

     Author: Jose Casadiego
     Date:   May 2017
    '''

    N,M = X.shape
    Expansion = np.zeros((K+1, M, N))

    if TYPE == 'polynomial':
        for n in xrange(N):
            for k in xrange(K):
                Expansion[k,:,n] = X[n,:]**k

    elif TYPE == 'polynomial_diff':
        Xi =np.zeros((N,M))
        for m in xrange(M):
            Xi[:,m] = X[:,m]-X[NODE,m]

        for n in xrange(N):
            for k in xrange(K):
                Expansion[k,:,n] = Xi[n,:]**k

    elif TYPE == 'fourier':
        Expansion=np.zeros((2*(K), M,N))
        for n in xrange(N):
            t = 0
            for k in xrange(K):
                Expansion[k+t,:,n] = np.sin(k*X[n,:])
                Expansion[k+t+1,:,n] = np.cos(k*X[n,:])
                t += 1

    elif TYPE == 'fourier_diff':
        Expansion=np.zeros((2*(K),M,N))
        Xi = np.zeros((N,M))
        for m in xrange(M):
            Xi[:,m] = X[:,m] - X[NODE, m]

        for n in xrange(N):
            t = 0
            for k in xrange(K):
                Expansion[k+t,:,n] = np.sin(k*Xi[n,:])
                Expansion[k+t+1,:,n] = np.cos(k*Xi[n,:])
                t += 1

    elif TYPE == 'power_series':
        Expansion = np.zeros(((K)*(K), M, N))
        for n in xrange(N):
            for k1 in xrange(K):
                for k2 in xrange(K):
                    for m in xrange(M):
                        Expansion[((K)*k1)+k2, m, n] = (X[NODE,m]**k1)*(X[n,m]**k2)

    elif TYPE == 'RBF':
        Expansion = np.zeros((K,M,N))
        for n in xrange(N):
            A = np.vstack((X[n,:], X[NODE,:]))
            for m1 in xrange(K):
                for m2 in xrange(M):
                    Expansion[m1,m2, n] = np.sqrt(2.0+np.linalg.norm(A[:,m1]-A[:,m2],2)**2)

    return(Expansion)


if __name__ == "__main__":
    import unittest
    from numpy.testing import assert_equal

    class TestExpansion(unittest.TestCase):
        def test_polynomial(self):
            X = np.reshape(np.arange(0,6,1), (3,2))
            K = 1
            NODE = 0
            E = basis_expansion(X,K,'polynomial',NODE)
            T = np.array([[[1.,1.,1.],
                          [1.,1.,1.]],
                        [[0., 2., 4.],
                         [1.,3.,5.,]]])
            assert_equal(E, T)

        def test_polynomial_diff(self):
            X = np.reshape(np.arange(0,6,1), (3,2))
            K = 1
            NODE = 0
            E = basis_expansion(X,K,'polynomial_diff',NODE)
            T = np.array([[[1.,1.,1.],
                          [1.,1.,1.]],
                        [[0., 2., 4.],
                         [0.,2.,4.,]]])
            assert_equal(E, T)

        def test_fourier(self):
            X = np.array([[(np.pi)/2., (np.pi)],
                          [(np.pi) , (np.pi)/2.]])
            K = 1
            NODE = 0
            T1 = np.array([[0., 0.],
                         [1., 1.],
                         [1., 0.],
                         [0., -1.]])

            T2 = np.array([[0., 0.],
                         [1., 1.],
                         [0., 1.],
                         [-1.,0.]])
            T = np.zeros((4,2,2))
            T[:,:,0] = T1
            T[:,:,1] = T2

            E = basis_expansion(X,K,'fourier',NODE)

            np.testing.assert_almost_equal(E, T)

        def test_fourier_diff(self):
            X = np.array([[(np.pi)/2., (np.pi)],
                          [(np.pi) , (np.pi)/2.]])
            K = 1
            NODE = 0
            E = basis_expansion(X,K,'fourier_diff',NODE)

            T1 = np.array([[0., 0.],
                         [1., 1.],
                         [0., 0.],
                         [1., 1.]])

            T2 = np.array([[0.,0.],
                           [1., 1.],
                           [1., -1.],
                           [0., 0.]])

            T = np.zeros((4,2,2))
            T[:,:,0] = T1
            T[:,:,1] = T2

            np.testing.assert_almost_equal(E, T)

        def test_power_series(self):
            X = np.reshape(np.arange(0,4,1), (2,2))
            K = 1
            NODE = 0
            E = basis_expansion(X,K,'power_series',NODE)

            T1 = np.array([[1., 1.],
                           [0., 1.],
                           [0., 1.],
                           [0., 1.]])
            T2 = np.array([[1., 1.],
                           [2., 3.],
                           [0., 1.],
                           [0., 3.]])
            T = np.zeros((4,2,2))
            T[:,:,0] = T1
            T[:,:,1] = T2

            np.testing.assert_almost_equal(E, T)

    unittest.main()
