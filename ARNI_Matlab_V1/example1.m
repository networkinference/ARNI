% example1.m generates different time series of networks of dynamical
% systems starting from different intial conditions and reconstructs the
% connectivity for a selected unit. Increasing the number of different time
% series leads to better results.
% 
% Parameters
% ------------------
% MODEL: Dynamical model on network units. Currently, only kuramoto1,
%        kuramoto2, michaelis_menten and roessler are supported. For
%        detailed information about the models, please check methods
%        section in the main manuscript.
% N:     Network size.
% NI:    Number of incoming connections per unit.
% S:     Number of different time series.
% M:     Number of time points per time series.
% NODE:  Unit upon the reconstruction takes place.
% BASIS: Type of basis employed. Currently, polynomial, polynomial_diff,
%        power_series, fourier, fourier_diff and RBF are supported. For
%        more detailed information, please see 'Functions/basis_expansion.m'
%        and Table I in the main manuscript.
% ORDER: Number of bases in the expansion.
%
% Input type
% ------------------
% MODEL: string
% N:     integer
% NI:    integer (NI<N)
% S:     integer
% M:     integer
% NODE:  integer
% BASIS: string
% ORDER: integer
%
% Output
% ------------------
% Figure showing: (1) evolution of fitting costs versus the number of inferred
% interactions with actual inferred interactions; and, (2) Receiver-Operating
% -Characteristic Curve.
%
% Accompanying material to "Model-free inference of direct interactions 
% from nonlinear collective dynamics".
%
% Author: Jose Casadiego
% Date:   May 2017

close all;
addpath('Models/','Functions/')

MODEL={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
BASIS={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};

N=25;
NI=4;
S=50;
M=5;
ORDER=6;

simulate(MODEL{1},N,NI,S,M);

NODE=21;

[list,cost,FPR,TPR,AUC]=reconstruct(MODEL{1},NODE,BASIS{2},ORDER);

figure('Name',['Reconstruction for unit ',num2str(NODE)]);
ax1=subplot(2,1,1);
plot(cost,'-o','LineWidth',2.5,'Color',[0,0.7,0.9],'MarkerFaceColor',[0,0.7,0.9]);
title({'Evolution of Fitting Costs';'(Actual connections shown above)'});
xlabel('# Inferred Interactions');
ylabel('Cost');

for i=1:length(cost)
    text(1.01*i,1.03*cost(i),num2str(list(i)))
end

ax2=subplot(2,1,2);
plot(FPR,TPR,'LineWidth',2.5,'Color',[0,0.7,0.9])
title('Receiver-Operating-Characteristic Curve');
xlabel('False Positives Rate');
ylabel('True Positives Rate');
text(0.4,0.5,['AUC score=',num2str(AUC)])
