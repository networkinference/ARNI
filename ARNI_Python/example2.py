#!/usr/bin/env python

from ARNIpy import simulate, reconstruct
import matplotlib.pyplot as plt

'''
 example2.m generates different time series for two different dynamical
 systems under different types of coupling functions, ($$h(x_j)$$ and $$h(x_i,x_j)$$)
 starting from different intial conditions and reconstructs the connectivity
 for a selected unit. Increasing the number of different time series leads
 to better results. Bivariate coupling functions are only correctly
 represented by bivariate coupling functions. Analogously, univariate
 functions are only correctly represented by univariate basis functions.

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
 interactions using different bases for models kuramoto2 and
 michaelis_menten.

 Accompanying material to "Model-free inference of direct interactions
 from nonlinear collective dynamics".

 Author: Jose Casadiego
 Date:   May 2017
'''

MODEL=('kuramoto1','kuramoto2','michaelis_menten','roessler')
BASIS=('polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF')
NAMES=('Polynomial','Polynomial Diff','Fourier','Fourier Diff','Power Series','Radial Basis Function')

N=25
NI=4
S=30
M=10
NODE=15
ORDER=6

simulate(MODEL[1],N,NI,S,M)

#This may take several minutes
f, axes = plt.subplots(2,3)
axes = axes.ravel()
st = f.suptitle('Reconstruction with different basis for ' + MODEL[1])
for i in range(len(BASIS)):
    llist, cost, FPR, TPR, AUC = reconstruct(MODEL[1], NODE, BASIS[i], ORDER)
    axes[i].plot(range(1,len(cost)+1), cost)
    axes[i].scatter(range(1,len(cost)+1), cost)
    axes[i].set_title('Fitting Costs: ' + NAMES[i] + '\nAUC=%.3f' % AUC,
    fontsize=8)
    axes[i].set_xlabel('# Inferred Interactions')
    axes[i].set_ylabel('Cost')

f.tight_layout()
st.set_y(0.95)
f.subplots_adjust(top=0.80)

#f.savefig('./Graphs/example2_kuramoto2.pdf')

simulate(MODEL[2],N,NI,S,M)

f, axes = plt.subplots(2,3)
axes = axes.ravel()
st = f.suptitle('Reconstruction with different basis for ' + MODEL[2])
for i in range(len(BASIS)):
    llist, cost, FPR, TPR, AUC = reconstruct(MODEL[2], NODE, BASIS[i], ORDER)
    axes[i].plot(range(1,len(cost)+1), cost)
    axes[i].scatter(range(1,len(cost)+1), cost)
    axes[i].set_title('Fitting Costs: ' + NAMES[i] + '\nAUC=%.3f' % AUC,
    fontsize=8)
    axes[i].set_xlabel('# Inferred Interactions')
    axes[i].set_ylabel('Cost')

f.tight_layout()
st.set_y(0.95)
f.subplots_adjust(top=0.80)

#f.savefig('./Graphs/example2_Mich_menten.pdf')
plt.show()
