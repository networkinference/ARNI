% example4.m compares the quality of reconstruction between short and long 
% time series with poor temporal resolution on networks of coupled
% chaotic roessler systems. Several short time series are preferable over
% long time series for this type of systems.
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
% Figures showing the evolution of fitting costs versus the number of inferred
% interactions using different number of bases on kuramoto2 models.
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
NODE=randperm(25,4);

N=25;
NI=4;
S=50;
M=5;
ORDER=6;

simulate(MODEL{4},N,NI,S,M);

figure('Name','Reconstruction of oscillators using short time series');

% This may take several minutes
t=1;
for node=NODE
    [list,cost,FPR,TPR,AUC]=reconstruct(MODEL{4},node,BASIS{1},ORDER);
    subplot(2,2,t);
    plot(FPR,TPR,'LineWidth',2.5,'Color',[0,0.7,0.9]);
    title(['ROC Curve-Unit: ',num2str(node)]);
    xlabel('FPR');
    ylabel('TPR');
    text(0.4,0.5,['AUC score=',num2str(AUC)])
    t=t+1;
end
set(gcf,'Position',[0 0 1000 500])

S=5;
M=50;

simulate(MODEL{4},N,NI,S,M);

figure('Name','Reconstruction of oscillators using long time series');

% This may take several minutes
t=1;
for node=NODE
    [list,cost,FPR,TPR,AUC]=reconstruct(MODEL{4},node,BASIS{1},ORDER);
    subplot(2,2,t);
    plot(FPR,TPR,'LineWidth',2.5,'Color',[0,0.7,0.9]);
    title(['ROC Curve-Unit: ',num2str(node)]);
    xlabel('FPR');
    ylabel('TPR');
    text(0.4,0.5,['AUC score=',num2str(AUC)])
    t=t+1;
end
set(gcf,'Position',[0 0 1000 500])