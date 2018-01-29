#!/usr/bin/env python

from ARNIpy import simulate, reconstruct
import matplotlib.pyplot as plt


''' example1.py generates different time series of networks of dynamical
 systems starting from different intial conditions and reconstructs the
 connectivity for a selected unit. Increasing the number of different time
 series leads to better results.

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
 NODE:  Unit upon the reconstruction takes place. Zero indexed
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
 Figure showing: (1) evolution of fitting costs versus the number of inferred
 interactions with actual inferred interactions; and, (2) Receiver-Operating
 -Characteristic Curve.

 Accompanying material to "Model-free inference of direct interactions
 from nonlinear collective dynamics".

 Author: Jose Casadiego
 Date:   May 2017
'''


MODEL=['kuramoto1', 'kuramoto2', 'michaelis_menten', 'roessler']
BASIS = ['polynomial', 'polynomial_diff', 'fourier', 'fourier_diff', 'power_series', 'RBF']

N=20
NI=4
S=10
M=10
ORDER=10

simulate(MODEL[0], N, NI, S, M)

NODE=0

llist, cost, FPR, TPR, AUC = reconstruct(MODEL[0], NODE, BASIS[5], ORDER)

f, (ax1,ax2) = plt.subplots(2,1)
st = f.suptitle('Reconstruction for unit ' + str(NODE))
ax1.plot(range(1,len(cost)+1), cost)
ax1.scatter(range(1,len(cost)+1), cost)
ax1.set_title('Evolution of fitting costs\n(Actual connections shown above)')
ax1.set_xlabel('# Inferred Interactions')
ax1.set_ylabel('Cost')

for i,c in enumerate(cost):
    ax1.text(1.01*(i+1),1.0*c, str(llist[i]))

ax2.plot(FPR, TPR)
ax2.set_title('Receiver-Operating-Characteristic Curve')
ax2.set_xlabel('False Positives Rate')
ax2.set_ylabel('True Positives Rate')
ax2.set_ylim((-0.1,1.1))
ax2.text(0.7, 0.1, 'AUC score=%.3f' % AUC)

f.tight_layout()
#
st.set_y(0.95)
f.subplots_adjust(top=0.8)
#f.savefig('../graphs/kuramoto1_fourier_diff.pdf')
plt.show()
