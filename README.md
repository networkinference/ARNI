# ARNI -- Algorithm for Revealing Network Interactions
Example codes in Matlab for simulating and reconstructing network dynamical systems.

## 1 Codes for simulating and reconstructing networks

Along this Supplementary Software, we provide a set of 
codes (including examples) for simulating time series and 
reconstructing the connectivity of different models of network 
dynamical systems. The codes we provide include:

ARNI_NatComm

├──Data

├──────connectivity.dat
   
├──────data.dat

├──────frequencies.dat

├──────ts_param.dat

├──Functions

├──────basis_expansion.m

├──────reconstruct.m

├──────simulate.m

├──────topology.m 

├──Html_Output

├──Models

├──────kuramoto1.m    
   
├──────kuramoto2.m 
   
├──────michaelis_menten.m   
   
├──────roessler.m

├──example1.m 

├──example2.m 

├──example3.m 

├──example4.m 

### 1.1 Data

The folder Data contains all the necessary information for 
reconstructing and evaluating the quality of reconstruction of 
networks from time series. Specifically, connectivity.dat 
contains the connectivity matrix of the network model simulated. 
The file data.dat contains all the simulated time series for a 
specific simulated model. In case of simulating phase-coupled 
oscillators, frequencies.dat contains the intrinsic frequency of 
each oscillator. Finally, ts_param.dat indicates how many time 
series and for how long such were simulated.

### 1.2 Functions

#### 1.2.1 basis_expansion.m 

basis_expansion.m generates a multidimensional array of basis 
expansions evaluated on all points of a multivariate time series.

##### Output

Multidimensional array containing the evalation of all basis 
functions for all time points and all possible incoming 
connections.

#### 1.2.2 reconstruct.m

reconstruct.m returns a ranked list of the inferred incoming 
connections.

##### Output

list: Sequence of inferred interactions (node identities) in the 
order such were detected.
cost: Fitting cost for all inferred interactions in the order 
such were detected.
FPR: False positives rate for the reconstruction.
TPR: True positives rate for the reconstruction.
AUC: Quality of reconstruction measured in AUC scores.

#### 1.2.3 simulate.m

simulate.m generates time series of networks of dynamical systems 
for several different intial conditions.

##### Output

'Data/data.dat': File containing all simulated time series in a 
concatenaded form.
'Data/ts_param.dat': File containing time series parameters, i.e. 
number and length of time series.

#### 1.2.4 topology.m

topology.m generates connectivity matrices for network 
simulation.

##### Output

'Data/connectivity.dat': File containing a weighted adjacency 
matrix.

We refer the reader to the headers of each function (and the 
examples described below) for more details about the proper usage 
of these functions.

## 1.3 Html_Output

This folder contains the examples and their outputs in HTML 
format for simple access.

## 1.4 Models

This folder contains the files for simulating random ensembles of gene regulatory circuits, phase-coupled oscillators with one and two fourier modes and networks of Rössler oscillators. 

## 1.5 Examples

example1.m generates different time series of networks of 
dynamical systems starting from different intial conditions and 
reconstructs the connectivity for a selected unit. Increasing the 
number of different time series leads to better results.

##### Output 

Figure showing: (1) evolution of fitting costs versus the number 
of inferred interactions with actual inferred interactions; and, 
(2) Receiver-Operating-Characteristic Curve.


example2.m generates different time series for two different 
dynamical systems under different types of coupling functions, h(x_{j})
 and h(x_{i},x_{j}), starting from different intial conditions 
and reconstructs the connectivity for a selected unit. Increasing 
the number of different time series leads to better results. 
Bivariate coupling functions are only correctly represented by 
bivariate coupling functions. Analogously, univariate functions 
are only correctly represented by univariate basis functions.

##### Output

Figures showing the evolution of fitting costs versus the number 
of inferred interactions using different bases for models 
kuramoto2 and michaelis_menten.


example3.m generates different time series for kuramoto2 systems 
and reconstructs them under radial basis functions of different 
orders. Greater orders (number of employed basis functions) lead 
to better results.

##### Output

Figures showing the evolution of fitting costs versus the number 
of inferred interactions using different number of bases.

example4.m compares the quality of reconstruction between short 
and long time series with poor temporal resolution on networks of 
coupled chaotic roessler systems. Several short time series are 
preferable over long time series for this type of systems.

##### Output

Figures showing the evolution of fitting costs versus the number 
of inferred interactions using different number of bases on 
kuramoto2 models.


