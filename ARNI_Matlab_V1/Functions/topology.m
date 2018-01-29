function topology(N,TYPE,DIRECTED,NI)
% topology(N,TYPE,DIRECTED,NI) generates connectivity matrices for network 
% simulation.
% 
% Parameters
% ------------------
% N:        Network size.
% TYPE:     Type of network. Currently, only homogeneous (random network
%           with fixed number of incoming connections) and regular are
%           supported.
% DIRECTED: Network (un)directionality, i.e. directed or undirected.
% NI:       Number of incoming connections per unit.
%
% Input type
% ------------------
% N:        integer
% TYPE:     string
% DIRECTED: string
% NI:       integer
%
% Output
% ------------------
% 'Data/connectivity.dat': File containing a weighted adjacency matrix.
%
% Example
% ------------------
% topology(20,'homogeneous','undirected',5) generates an undirected 
% connectivity matrix of 20x20, where each unit has 5 randomly-selected 
% incoming connections.
%
% Author: Jose Casadiego
% Date:   May 2017

types={'homogeneous','regular'};
directness={'directed','undirected'};
if any(ismember(types,TYPE))==0
    disp('ERROR: TYPE must be a valid string: homogeneous, regular');
    
elseif any(ismember(directness,DIRECTED))==0
    disp('ERROR: DIRECTED must be a valid string: directed, undirected');
    
else
    J=zeros(N,N);           %coupling matrix
    
    switch TYPE
            
        case 'homogeneous'  %homogeneous topology with NI connections per unit
            
            for i=1:N
                f=randperm(N);
                f(f==i)=[];
                for j=1:NI
                    J(i,f(j))=(0.5+(1-0.5)*rand(1,1))/NI;
                end
            end
            
        case 'regular'      %regular structure with NI connections per unit
            
            f=randperm(N);
            f(f==1)=[];
            
            a=zeros(1,N);
            for i=1:NI
                a(1,f(i))=1;
            end
            
            for n=1:N
                J(n,:)=a;
                b=a(1,end);
                a=a(1,1:end-1);
                a=horzcat(b,a);
            end
            
            R=(0.5+(2-0.5)*rand(N))/NI;
            
            J=J.*R;

    end
    
    switch DIRECTED
        
        case 'directed'
            
            dlmwrite('Data/connectivity.dat', J, 'delimiter', '\t', 'precision', 4);
            
        case 'undirected'
            
            for i=1:N
                for j=1:N
                    if J(i,j)~=0 && J(j,i)==0
                        J(j,i)=J(i,j);
                    end
                end
            end
            
            dlmwrite('Data/connectivity.dat', J, 'delimiter', '\t', 'precision', 4);
            
    end
    
end
end