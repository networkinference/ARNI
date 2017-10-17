function [list,cost,FPR,TPR,AUC]=reconstruct(MODEL,NODE,BASIS,ORDER)
% reconstruct(MODEL,NODE,BASIS,ORDER) returns a ranked list of the
% inferred incoming connections.

% Parameters
% ------------------
% MODEL: Dynamic model employed. This is only used to specify whether the
%        time series come from 1D systems like kuramoto1 or 3D systems like
%        roessler. Thus, it is not used during the actual reconstruction.
% NODE:  Unit upon the reconstruction takes place.
% BASIS: Type of basis employed. Currently, polynomial, polynomial_diff,
%        power_series, fourier, fourier_diff and RBF are supported. For
%        more detailed information, please see 'Functions/basis_expansion.m'
%        and Table I in the main manuscript.
% ORDER: Number of basis in the expansion.
%
% Input type
% ------------------
% MODEL: string
% NODE:  integer
% BASIS: string
% ORDER: integer
%
% Output
% ------------------
% list: Sequence of inferred interactions in the order such were detected.
% cost: Fitting cost for all inferred interactions in the order such were
%       detected.
% FPR:  False positives rate for the reconstruction.
% TPR:  True positives rate for the reconstruction.
% AUC:  Quality of reconstruction measured in AUC scores.
%
% Example
% ------------------
% reconstruct('michaelis_menten',10,'polynomial',6) reconstructs the
% connectivity of unit 10 using polynomials up to power 6 as basis functions.
%
% Accompanying material to "Model-free inference of direct interactions 
% from nonlinear collective dynamics".
%
% Author: Jose Casadiego
% Date:   May 2017

% Stopping criterium: decrease it to recover longer list of possible links
th=0.0001;

models={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
bases={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};

if any(ismember(models,MODEL))==0
    disp('ERROR: MODEL must be a valid string: kuramoto1, kuramoto2, michaelis_menten, roessler');
elseif any(ismember(bases,BASIS))==0
    disp('ERROR: BASIS must be a valid string: polynomial, polynomial_diff, fouries, fourier_diff, power_series, RBF');
    
else
    disp('Initiating reconstruction...');
    disp('Reading data...');
    data = dlmread('Data/data.dat');
    connectivity=dlmread('Data/connectivity.dat');
    ts_param=dlmread('Data/ts_param.dat');
    data=data';
    
    S=ts_param(1,1);
    M=ts_param(1,2);
    [N,~]=size(data);
    
    x=data;
    
    % Estimating time derivatives and constructing input matrices
    disp('Estimating time derivatives and constructing input matrices...');
    Xtemp=[];
    DX=[];
    for s=1:S
        Ytemp=zeros(N,M);
        DY=zeros(N,M);
        for n=1:N
            for m=1:M
                Ytemp(n,m)=(x(n,m+(s-1)*(M+1))+x(n,m+1+(s-1)*(M+1)))*0.5;
                DY(n,m)=(-x(n,m+(s-1)*(M+1))+x(n,m+1+(s-1)*(M+1)))*1/(1.0);
            end
        end
        Xtemp=[Xtemp,Ytemp];
        DX=[DX,DY];
    end
    
    switch MODEL
        
        case 'roessler'
            
            % Construction of connectivity matrix for Roessler oscillators
            % including y and z variables
            Ns=ceil(N/3);
            connectivity2=zeros(Ns,N);
            
            for i=1:Ns
                for j=1:Ns
                    connectivity2(i,3*(j-1)+1)=connectivity(i,j);
                    if i==j
                        connectivity2(i,3*(j-1)+2)=1;
                        connectivity2(i,3*(j-1)+3)=1;
                    end
                end
            end
            
            X=Xtemp;
            
            % Beginning of reconstruction algorithm
            disp('Performing ARNI...');
            Y=basis_expansion(X,ORDER,BASIS,NODE);
            nolist=1:N;
            list=[];
            cost=[];
            b=1;
            vec=zeros(1,N);
            while (~isempty(nolist)) && (b==1);
                % Composition of inferred subspaces
                Z=[];
                for n=1:length(list)
                    Z=[Z;Y(:,:,list(n))];
                end
                
                % Projection on remaining composite spaces
                P=zeros(length(nolist),2);
                cost_err=zeros(length(nolist),1);
                for n=1:length(nolist)
                    % Composition of a possible space
                    R=[Z;Y(:,:,nolist(n))];
                    % Error of projection on possible composite space
                    P(n,1)=std(DX(3*(NODE-1)+1,:)-DX(3*(NODE-1)+1,:)*pinv(R)*R);
                    P(n,2)=nolist(n);
                    % Fitting cost of possible composite space
                    cost_err(n,1)=1/M *norm(DX(3*(NODE-1)+1,:)-DX(3*(NODE-1)+1,:)*pinv(R)*R);
                    R=[];
                end
                
                if std(P(:,1))<th
                    b=0;
                    break
                    
                else
                    % Selection of composite space which minimizes
                    % projection error
                    [MIN,block]=min(P(:,1));
                    list=[list,P(block,2)];
                    nolist(nolist==P(block,2))=[];
                    vec(1,P(block,2))=MIN;
                    cost=[cost,cost_err(block,1)];
                end
                
            end
            % End of reconstruction algorithm
            
            adjacency=connectivity2;
            adjacency(adjacency~=0)=1;
            
            % Evaluation of results via AUC score
            [FPR,TPR,~,AUC]=perfcurve(abs(adjacency(NODE,:)),abs(vec),1);
            
            disp('Reconstruction has finished!');
            disp('Quality of reconstruction:');
            disp(AUC)
                 
        otherwise
            
            if (strcmp(MODEL,'kuramoto1')) || (strcmp(MODEL,'kuramoto2'))
                % Transforming data coming from phase oscillators
                X=mod(Xtemp,2*pi);
            else
                X=Xtemp;
            end
            
            % Beginning of reconstruction algorithm
            disp('Performing ARNI...');
            Y=basis_expansion(X,ORDER,BASIS,NODE);
            nolist=1:N;
            list=[];
            cost=[];
            b=1;
            vec=zeros(1,N);
            while (~isempty(nolist)) && (b==1);
                % Composition of inferred subspaces
                Z=[];
                for n=1:length(list)
                    Z=[Z;Y(:,:,list(n))];
                end
                
                % Projection on remaining composite spaces
                P=zeros(length(nolist),2);
                cost_err=zeros(length(nolist),1);
                for n=1:length(nolist)
                    % Composition of a possible space
                    R=[Z;Y(:,:,nolist(n))];
                    % Error of projection on possible composite space
                    P(n,1)=std(DX(NODE,:)-DX(NODE,:)*pinv(R)*R);
                    P(n,2)=nolist(n);
                    % Fitting cost of possible composite space
                    cost_err(n,1)=1/M *norm(DX(NODE,:)-DX(NODE,:)*pinv(R)*R);
                    R=[];
                end
                
                if std(P(:,1))<th
                    b=0;
                    break
                else
                    % Selection of composite space which minimizes
                    % projection error
                    [MIN,block]=min(P(:,1));
                    list=[list,P(block,2)];
                    nolist(nolist==P(block,2))=[];
                    vec(1,P(block,2))=MIN;
                    cost=[cost,cost_err(block,1)];
                end
                
            end
            % End of reconstruction algorithm
            
            adjacency=connectivity;
            adjacency(adjacency~=0)=1;
            
            % Adding degradation rate to true adjacency matrix of
            % Michaelis-Menten systems
            if strcmp(MODEL,'michaelis_menten')
                adjacency(1:N+1:N*N) = 1;
            end
            
            % Evaluation of results via AUC score
            [FPR,TPR,~,AUC]=perfcurve(abs(adjacency(NODE,:)),abs(vec),1);
            disp('Reconstruction has finished!');
            disp('Quality of reconstruction:');
            disp(AUC);
            
    end
end