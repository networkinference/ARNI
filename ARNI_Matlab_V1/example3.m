% example3.m generates different time series for kuramoto2 systems and 
% reconstructs them under radial basis functions of different orders.
% Greater orders (number of employed basis functions) lead to better results.
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
% interactions using different number of bases.
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
S=30;
M=10;
NODE=15;

simulate(MODEL{2},25,4,30,10);

ORDER=[5,10,15,20,25,30];
figure('Name','Reconstruction using different number of RBF');

% This may take several minutes
t=1;
for k=ORDER
    [list,cost,FPR,TPR,AUC]=reconstruct(MODEL{2},NODE,BASIS{6},k);
    subplot(2,3,t);
    plot(cost,'-o','LineWidth',2.5,'Color',[0,0.7,0.9],'MarkerFaceColor',[0,0.7,0.9]);
    title({['Fitting Costs: ',num2str(k),' RBF'],['AUC=',num2str(AUC)]});
    xlabel('# Inferred Interactions');
    ylabel('Cost');
    t=t+1;
end
set(gcf,'Position',[0 0 1000 500])