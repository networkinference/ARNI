#!/usr/bin/env python

from ARNIpy import simulate, reconstruct
import matplotlib.pyplot as plt
'''
 example3.m generates different time series for kuramoto2 systems and
 reconstructs them under radial basis functions of different orders.
 Greater orders (number of employed basis functions) lead to better results.

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
 interactions using different number of bases.

 Accompanying material to "Model-free inference of direct interactions
 from nonlinear collective dynamics".

 Author: Jose Casadiego
 Date:   May 2017
'''

MODEL=('kuramoto1','kuramoto2','michaelis_menten','roessler')
BASIS=('polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF')

N=25
NI=4
S=30
M=10
NODE=15

simulate(MODEL[1],N,NI,S,M)
ORDER=(5,10,15,20,25,30)
f, axes = plt.subplots(2,3)
st = f.suptitle('Reconstruction with different number of RBF')

axes = axes.ravel()
#this may take several minutes
for t,k in enumerate(ORDER):
    llist, cost, FPR, TPR, AUC = reconstruct(MODEL[0], NODE, BASIS[5], k)
    axes[t].plot(range(1,len(cost)+1), cost)
    axes[t].scatter(range(1,len(cost)+1), cost)
    axes[t].set_title('Fitting costs: ' + str(k) + ' RBF\nAUC=%.3f' % AUC, fontsize=8)
    axes[t].set_xlabel('# Inferred Interactions')
    axes[t].set_ylabel('Cost')

f.tight_layout()
st.set_y(0.95)
f.subplots_adjust(top=0.80)
plt.show()

#f.savefig('./Graphs/example3.pdf')
