#!/usr/bin/env python
from math import pi
import numpy as np
from sklearn.metrics import roc_curve, auc
import sys

from basis_expansion import basis_expansion

def reconstruct(MODEL, NODE, BASIS, ORDER):
    '''
    reconstruct(MODEL, NODE, BASIS, ORDER) returns a ranked list of the inferred
    incoming connections

     Parameters
     ------------------
     MODEL: Dynamic model employed. This is only used to specify whether the
            time series come from 1D systems like kuramoto1 or 3D systems like
            roessler. Thus, it is not used during the actual reconstruction.
     NODE:  Unit upon the reconstruction takes place. Zero indexed
     BASIS: Type of basis employed. Currently, polynomial, polynomial_diff,
            power_series, fourier, fourier_diff and RBF are supported. For
            more detailed information, please see 'Functions/basis_expansion.m'
            and Table I in the main manuscript.
     ORDER: Number of basis in the expansion.

     Input type
     ------------------
     MODEL: string
     NODE:  integer
     BASIS: string
     ORDER: integer

     Output
     ------------------
     list: Sequence of inferred interactions in the order such were detected.
     cost: Fitting cost for all inferred interactions in the order such were
           detected.
     FPR:  False positives rate for the reconstruction.
     TPR:  True positives rate for the reconstruction.
     AUC:  Quality of reconstruction measured in AUC scores.

     Example
     ------------------
     reconstruct('michaelis_menten',10,'polynomial',6) reconstructs the
     connectivity of unit 10 using polynomials up to power 6 as basis functions.

     Accompanying material to "Model-free inference of direct interactions
     from nonlinear collective dynamics".

     Author: Jose Casadiego
     Date:   May 2017
     '''

    #Stopping criterium: decrease it to recover longer list of possible links
    th=0.0001

    models=['kuramoto1', 'kuramoto2', 'michaelis_menten', 'roessler']
    bases=['polynomial', 'polynomial_diff', 'fourier', 'fourier_diff', 'power_series', 'RBF']

    if (MODEL not in models):
        sys.exit('ERROR: MODEL must be a valid string: kuramoto1')

    elif (BASIS not in bases):
        sys.exit('ERROR: BASIS must be a valid string: polynomial')

    else:
        print('Initiating reconstruction...')
        print('Reading data...')
        data = np.loadtxt('Data/data.dat', delimiter='\t')
        connectivity = np.loadtxt('Data/connectivity.dat', delimiter='\t')
        ts_param=np.loadtxt('Data/ts_param.dat', delimiter='\t')
        data=data.transpose()

        S = int(ts_param[0])
        M = int(ts_param[1])
        N = int(data.shape[0])

        x = data

        # Estimating time derivatives and constructing input matrices
        print('Estimating time derivatives and constructing input matrices...')
        Xtemp = np.array([])
        DX = np.array([])
        for s in xrange(S):
            m_start = M*s
            m_end = M*s + (M-1)
            x0 = x[:,m_start:m_end]
            x1 = x[:,m_start+1:m_end+1]

            Ytemp = (x0 + x1) * 0.5
            DY = x1 - x0

            Xtemp = np.hstack((Xtemp, Ytemp)) if Xtemp.size else Ytemp
            DX = np.hstack((DX, DY)) if DX.size else DY

        if MODEL == 'roessler':
            # Construction of connectivity matrix for Roessler oscillators
            # including y and z variables
            Ns = int(np.ceil(N/3.))
            connectivity2 = np.zeros((Ns,N))

            # for each node, x1 regulated by x2 & x3. Modify adjacency matrix
            # for comparison to predictions
            for i in xrange(Ns):
                for j in xrange(Ns):
                    connectivity2[i,(3*j)+0] = connectivity[i,j]
                    if i == j:
                        connectivity2[i,(3*j)+1] = 1
                        connectivity2[i, (3*j)+2] = 1

            X = Xtemp

            #beginning of reconstruction algorithm
            print('Performing ARNI...')
            Y = basis_expansion(X, ORDER, BASIS, NODE)
            nolist = range(N)
            llist = []
            cost = []
            b = 1
            vec = np.zeros(N,)
            while (nolist and (b==1)):
                #composition of inferred subspaces
                Z = np.array([])
                for n in xrange(len(llist)):
                    Z = np.vstack((Z, Y[:,:,llist[n]])) if Z.size else Y[:,:,llist[n]]

                # Projection on remaining composite subspaces
                P = np.zeros((len(nolist), 2))
                cost_err = np.zeros(len(nolist),)
                for n in xrange(len(nolist)):
                    #composition of a possible space
                    R = np.vstack((Z, Y[:,:,nolist[n]])) if Z.size else Y[:,:,nolist[n]]
                    #Error of projection on possible composite space
                    #  ( A.R=DX )
                    RI = np.linalg.pinv(R)
                    A = np.dot(DX[(3*NODE),:], RI)
                    DX_est = np.dot(A, R)
                    DIFF = DX[(3*NODE),:] - DX_est
                    P[n,0] = np.std(DIFF)
                    P[n,1] = int(nolist[n])
                    # Fitting cost of possible composite
                    cost_err[n] = (1/float(M)) * np.linalg.norm(DIFF)
                    R = np.array([])

                if np.std(P[:,0]) < th:
                    b = 0
                    break

                else:
                    # Selection of composite space which minimizes projection error
                    MIN = np.min(P[:,0]) #best score
                    block = np.argmin(P[:,0]) #node index of best
                    llist.append(int(P[block,1])) # add best node ID to llist
                    nolist.remove(int(P[block,1])) # remove best from candidate list
                    vec[int(P[block,1])] = MIN # used in ROC curve
                    cost.append(cost_err[block]) # record SS error

            # end of reconstruction algorithm
            if not llist:
                print('WARNING: no predicted regulators - check that NODE abundance varies in the data!')
                AUC = np.nan
                FPR = [np.nan]
                TPR = [np.nan]

            #evaluation of results via AUC score
            else:
                #load connectivity for comparison
                adjacency = connectivity2
                adjacency[adjacency != 0] = 1

                print('Quality of reconstruction:')

                if (np.sum(adjacency[NODE,:]) == 0):
                    print('WARNING: no true positive regulators!')
                    AUC = np.nan
                    FPR = [np.nan]
                    TPR = [np.nan]
                else:
                    # Evaluation of results via AUC score
                    FPR, TPR, _ = roc_curve(np.abs(adjacency[NODE,:]), np.abs(vec), 1)
                    AUC = auc(FPR, TPR)

                    FPR = np.insert(FPR,0,0.)
                    TPR = np.insert(TPR,0,0.)

                print(AUC)



        else: #if not Roessler
            if MODEL in ('kuramoto1', 'kuramoto2'):
                #Transforming data coming from phase oscillators
                X = np.mod(Xtemp, 2*pi)
            else:
                X = Xtemp

            # Beginning of reconstruction algorithm
            print('Performing ARNI...')
            # Y[basis, sample, node]
            Y = basis_expansion(X, ORDER, BASIS, NODE)
            nolist = range(N)
            llist = []
            cost = []
            b=1
            vec = np.zeros(N,)
            while (nolist and (b==1)):
                #composition of inferred subspaces
                Z = np.array([])
                for n in xrange(len(llist)):
                    Z = np.vstack((Z,Y[:,:,llist[n]])) if Z.size else Y[:,:,llist[n]]

                # projection on remaining composite spaces
                P = np.zeros((len(nolist),2))
                cost_err = np.zeros(len(nolist),)
                for n in xrange(len(nolist)):
                    #composition of a possible spaces
                    R = np.vstack((Z, Y[:,:,nolist[n]])) if Z.size else Y[:,:,nolist[n]]
                    #error of projection on possible composite space
                    # ( A.R=DX)
                    RI = np.linalg.pinv(R)
                    A = np.dot(DX[NODE,:], RI)
                    DX_est = np.dot(A, R)
                    DIFF = DX[NODE,:] - DX_est
                    P[n,0] = np.std(DIFF) # the uniformity of error
                    P[n,1] = int(nolist[n])
                    #Fitting cost of possible composite space
                    cost_err[n] =  (1/float(M)) * np.linalg.norm(DIFF)
                    R = np.array([])

                # break if all candidates equivalent
                if np.std(P[:,0]) < th:
                    b=0
                    break

                else:
                    #Selection of composite space which minimises projection error
                    MIN = np.min(P[:,0]) #best score
                    block = np.argmin(P[:,0]) #node index of best
                    llist.append(int(P[block,1])) # add best node ID to llist
                    nolist.remove(int(P[block,1])) # remove best from candidate list
                    vec[int(P[block,1])] = MIN # used in ROC curve
                    cost.append(cost_err[block]) # record SS Error
            print('Reconstruction has finished!')

            if not llist:
                print('WARNING: no predicted regulators - check that NODE abundance varies in the data!')
                AUC = np.nan
                FPR = [np.nan]
                TPR = [np.nan]

            else:
                #load connectivity for comparison
                adjacency = connectivity
                adjacency[adjacency != 0] = 1

                #adding degradation rate to true adjecency matric of Michaelis-menten systems
                if MODEL == 'michaelis_menten':
                    for i in xrange(N):
                        adjacency[i,i] = 1

                print('Quality of reconstruction:')

                if (np.sum(adjacency[NODE,:]) == 0):
                    print('WARNING: no true regulators!')
                    AUC = np.nan
                    FPR = [np.nan]
                    TPR = [np.nan]
                else:
                    # Evaluation of results via AUC score
                    FPR, TPR, _ = roc_curve(np.abs(adjacency[NODE,:]),
                    np.abs(vec), 1)
                    AUC = auc(FPR, TPR)
                    FPR = np.insert(FPR,0,0.)
                    TPR = np.insert(TPR,0,0.)

                print(AUC)

        return(llist, cost, FPR, TPR, AUC)

if __name__ =='__main__':
    llist, cost,_,_,_= reconstruct('roessler', 0,'polynomial', 2)
    print(llist)
    print(cost)
