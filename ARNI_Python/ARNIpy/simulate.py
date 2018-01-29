#!/usr/bin/env python
import numpy as np
from scipy.integrate import odeint
import subprocess
import sys

from models import kuramoto1,kuramoto2, michaelis_menten, roessler
from topology import topology

def simulate(MODEL, N, NI, S, M):
    '''
     simulate(MODEL,N,NI,S,M) generates time series of networks of dynamical
     systems for several different intial conditions.

     Parameters
     ------------------
     MODEL: Dynamical model on network units. Currently, only kuramoto1,
            kuramoto2, michaelis_menten and roessler are supported. For
            detailed information about the models, please check methods
            section in the main manuscript.
     N:     Network size.
     NI:    Number of incoming connections per unit.
     S:     Number of different time series.
     M:     Number of time points per time series.

     Input type
     ------------------
     MODEL: string
     N:     integer
     NI:    integer (NI<N)
     S:     integer
     M:     integer

     Output
     ------------------
     'Data/data.dat':     File containing all simulated time series in a
                          concatenaded form.
     'Data/ts_param.dat': File containing time series parameters, i.e. S and
                          M, for later extracting the different time series.

     Example
     ------------------
     simulate('kuramoto2',25,4,30,10) generates 30 time series of 10 time
     points each for a network of 25 oscillators defined by the model
     kuramoto2. Each oscillator has 4 incoming connections.

     Accompanying material to "Model-free inference of direct interactions
     from nonlinear collective dynamics".

     Author: Jose Casadiego
     Date:   May 2017
    '''

    #smpling rate of time series
    resolution=1

    models={'kuramoto1', 'kuramoto2', 'michaelis_menten', 'roessler'}
    if (MODEL not in models):
        sys.exit('ERROR: MODEL must be a valid string:kuramoto1, kuramoto2, michaelis_menten, roessler')
    else:
        cmd='rm -r Data/'
        process = subprocess.Popen(cmd.split())
        process.wait()

        cmd='mkdir Data/'
        process = subprocess.Popen(cmd.split())
        process.wait()

        print('Creating network structure...')
        topology(N, 'homogeneous', 'directed', NI) #
        print('Simulating time series...')
        Y = np.array([])
        if MODEL == 'kuramoto1':
            w = -2 + 4*np.random.uniform(low=0., high=1., size=(N,))
            np.savetxt('Data/frequencies.dat', w, fmt='%.4f',delimiter='\t')
            for s in xrange(S):
                init = -3.14 + (3.14+3.14) * np.random.uniform(0.,1.,size=(N,))
                tspan = np.arange(0,M,resolution)
                y = odeint(kuramoto1, init, tspan)
                Y = np.vstack((Y,y)) if Y.size else y

        elif MODEL == 'kuramoto2':
            w = -2 + (4) * np.random.uniform(0.,1.,size=(N,))
            np.savetxt('Data/frequencies.dat',w, fmt='%.4f', delimiter='\t')
            for s in xrange(S):
                init = -3.14 + (3.14+3.14)*np.random.uniform(0.,1.,size=(N,))
                tspan=np.arange(0,M,resolution)
                y = odeint(kuramoto2, init, tspan)
                Y = np.vstack((Y,y)) if Y.size else y

        elif MODEL == 'michaelis_menten':
            for s in xrange(S):
                init = 1+np.random.uniform(0.,1.,size=(N,))
                tspan=np.arange(0,M,resolution)
                y = odeint(michaelis_menten, init, tspan)
                Y = np.vstack((Y,y)) if Y.size else y

        elif MODEL == 'roessler':
            for s in xrange(S):
                init=-5 + (5+5)*np.random.uniform(0.,1., size=(3*N,))
                tspan = np.arange(0,M,resolution)
                y = odeint(roessler, init, tspan)
                Y = np.vstack((Y,y)) if Y.size else y

        ts_param=[S,M]
        np.savetxt('Data/data.dat', Y, fmt='%.4f', delimiter='\t')
        np.savetxt('Data/ts_param.dat', ts_param, fmt='%i',delimiter='\t')

        print('Simulation finished!')

if __name__ =='__main__':
    simulate('michaelis_menten', 5, 2, 3, 10)
