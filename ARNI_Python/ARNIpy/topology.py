#!/usr/bin/env python
import numpy as np
import sys

def topology(N,TYPE, DIRECTED, NI):
    '''
     topology(N,TYPE,DIRECTED,NI) generates connectivity matrices for network
     simulation.

     Parameters
     ------------------
     N:        Network size.
     TYPE:     Type of network. Currently, only homogeneous (random network
               with fixed number of incoming connections) and regular are
               supported.
     DIRECTED: Network (un)directionality, i.e. directed or undirected.
     NI:       Number of incoming connections per unit.

     Input type
     ------------------
     N:        integer
     TYPE:     string
     DIRECTED: string
     NI:       integer

     Output
     ------------------
     'Data/connectivity.dat': File containing a weighted adjacency matrix.

     Example
     ------------------
     topology(20,'homogeneous','undirected',5) generates an undirected
     connectivity matrix of 20x20, where each unit has 5 randomly-selected
     incoming connections.

     Author: Jose Casadiego
     Date:   May 2017
    '''

    types={'homogeneous', 'regular'}
    directness={'directed', 'undirected'}
    if (TYPE not in types):
        sys.exit('ERROR: TYPE must be a valid string: homogeneous, regular')
    elif (DIRECTED not in directness):
        sys.exit('ERROR: DIRECTED must be a valid string: directed, undirected')
    else:
            J = np.zeros((N,N)) # coupling matrices

            if TYPE == 'homogeneous': #homogeneous topology with NI connection per unit
                for i in xrange(N):
                    f = np.random.permutation(range(N))
                    f = [x for x in f if x != i]
                    for j in xrange(NI):
                        J[i,f[j]]=(0.5+(1-0.5)*np.random.uniform(0,1))/NI

            elif TYPE == 'regular': #regular structure with NI connections per unit

                f = np.random.permutation(range(1,N))

                a = [0]*N
                for i in xrange(NI):
                    a[f[i]] = 1

                for n in xrange(N):
                    J[n,:] = a
                    b = a[-1]
                    a = a[:-1]
                    a = [b] + a

                R = (0.5+(1-0.5)*np.random.uniform(0,1,size=(N,N)))/NI
                J = J * R

            if DIRECTED == 'directed':
                np.savetxt('Data/connectivity.dat', J, fmt='%.4f', delimiter='\t')

            elif DIRECTED == 'undirected':
                for i in xrange(N):
                    for j in xrange(N):
                        if (J[i,j]!=0 and J[j,i] == 0):
                            J[j,i]=J[i,j]

                np.savetxt('Data/connectivity.dat', J, fmt='%.4f', delimiter='\t')

if __name__ == '__main__':
    topology(10, 'regular', 'directed', 4)
