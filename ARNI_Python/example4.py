#!/usr/bin/env python

from ARNIpy import simulate, reconstruct
import matplotlib.pyplot as plt
import numpy as np


'''
 example4.m compares the quality of reconstruction between short and long
 time series with poor temporal resolution on networks of coupled
 chaotic roessler systems. Several short time series are preferable over
 long time series for this type of systems.

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
 NODE:  Unit upon the reconstruction takes place.
 BASIS: Type of basis employed. Currently, polynomial, polynomial_diff,
        power_series, fourier, fourier_diff and RBF are supported. For
        more detailed information, please see 'Functions/basis_expansion.m'
        and Table I in the main manuscript.
 ORDER: Number of bases in the expansion.

 Input type
 ------------------
 MODEL: string
 N:     integer
 NI:    integer (NI<N)
 S:     integer
 M:     integer
 NODE:  integer
 BASIS: string
 ORDER: integer

 Output
 ------------------
 Figures showing the evolution of fitting costs versus the number of inferred
 interactions using different number of bases on kuramoto2 models.

 Accompanying material to "Model-free inference of direct interactions
 from nonlinear collective dynamics".

 Author: Jose Casadiego
 Date:   May 2017
'''

MODEL=('kuramoto1','kuramoto2','michaelis_menten','roessler')
BASIS=('polynomial','polynomial_diff','fourier','fourier_diff','power_series',
        'RBF')

N=25
NI=4
S=50
M=5
ORDER=6
x = range(N)
np.random.shuffle(x)
NODE = x[:4]

simulate(MODEL[3],N,NI,S,M)

f1,axes1 = plt.subplots(2,2)
st1 = f1.suptitle('Reconstruction of oscillators using short time series')
axes1 = np.ravel(axes1)
#this may take several minutes
for t,node in enumerate(NODE):
    llist,cost, FPR, TPR, AUC = reconstruct(MODEL[3], node, BASIS[0], ORDER)
    axes1[t].plot(FPR,TPR)
    axes1[t].set_title('ROC-curve unit: ' + str(node), fontsize=8)
    axes1[t].set_xlabel('FPR')
    axes1[t].set_ylabel('TPR')
    axes1[t].text(0.5, 0.1, 'AUC score=%.3f' % AUC, fontsize=6)
f1.tight_layout()
st1.set_y(0.95)
f1.subplots_adjust(top=0.80)
#plt.savefig('./Graphs/example4_short_timeseries.pdf')
#plt.show()

S=5
M=50
simulate(MODEL[3],N,NI,S,M)

f2,axes2 = plt.subplots(2,2)
st2 = f2.suptitle('Reconstruction of oscillators using long time series')
axes2 = np.ravel(axes2)
#this may take several minutes
for t, node in enumerate(NODE):
    llist,cost,FPR,TPR,AUC = reconstruct(MODEL[3], node, BASIS[0], ORDER)
    axes2[t].plot(FPR,TPR)
    axes2[t].set_title('ROC-curve unit: ' + str(node), fontsize=8)
    axes2[t].set_xlabel('FPR')
    axes2[t].set_ylabel('TPR')
    axes2[t].text(0.5,0.1, 'AUC score=%.3f' % AUC, fontsize=6)
f2.tight_layout()
st2.set_y(0.95)
f2.subplots_adjust(top=0.80)
#plt.savefig('./Graphs/example4_long_timeseries.pdf')
plt.show()
