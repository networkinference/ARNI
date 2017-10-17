function Expansion = basis_expansion(X,K,TYPE,NODE)
% basis_expansion(X,K,TYPE,NODE) generates a multidimensional array of
% basis expansions evaluated on all points of a multivariate time series.
% 
% Parameters
% ------------------
% X:    Matrix containing N time series of M time points.
% K:    Maximum order of the basis expansion.
% TYPE: Type of basis function employed. In this file, we only
%       implement expansions up to pairwise interactions. Expansions
%       availables are: polynomial (x_{j}^{k}), polynomial_diff
%       ((x_{j}-x_{NODE})^{k}), power_series (x_{j}^{k1}*x_{i}^{k2}),
%       fourier (sin(k*x_{j}) and cos(k*x_{j})), fourier_diff 
%       (sin(k*(x_{j}-x_{i})) and cos(k*(x_{j}-x_{i}))) and RBF (a model
%       based on radial basis functions). These functions are shown in 
%       table I in the main manuscript. 
% NODE: Unit on which we are performing the reconstruction.
%
% Input type
% ------------------
% X:    double
% K:    integer
% TYPE: string
% NODE: integer
%
% Output
% ------------------
% Expansion: Multidimensional array of size [K+1,M,N] containing the
% evalation of all k=0,1,...,K basis functions for all M time points and 
% all N possible incoming connections. For power_series, (K*K+2) basis
% functions are employed, and for fourier(_diff), 2*(K+1) are employed. 
%
% Example
% ------------------
% basis_expansion(X,4,'power_series',5); generates a multidimensional array
% of size [18,M,N] containing the evaluation of the basis for all M time
% points and all N possible incoming connections.
%
% Accompanying material to "Model-free inference of direct interactions 
% from nonlinear collective dynamics".
%
% Author: Jose Casadiego
% Date:   May 2017

[N,M]=size(X);

Expansion=zeros(K+1,M,N);

switch TYPE
    
    case 'polynomial'
        
        for n=1:N
            for k=0:K
                Expansion(k+1,:,n)=X(n,:).^(k);
            end
        end
        
    case 'polynomial_diff'
        
        Xi=zeros(N,M);
        for m=1:M
            Xi(:,m)=X(:,m)-X(NODE,m);
        end
        
        for n=1:N
            for k=0:K
                Expansion(k+1,:,n)=Xi(n,:).^(k);
            end
        end
        
    case 'fourier'
        Expansion=zeros(2*(K+1),M,N);
        for n=1:N
            t=1;
            for k=0:K
                Expansion(k+t,:,n)=sin(k*X(n,:));
                Expansion(k+t+1,:,n)=cos(k*X(n,:));
                t=t+1;
            end
        end
    case 'fourier_diff'
        Expansion=zeros(2*(K+1),M,N);
        Xi=zeros(N,M);
        for m=1:M
            Xi(:,m)=X(:,m)-X(NODE,m);
        end
        
        for n=1:N
            t=1;
            for k=0:K
                Expansion(k+t,:,n)=sin(k*Xi(n,:));
                Expansion(k+t+1,:,n)=cos(k*Xi(n,:));
                t=t+1;
            end
        end
        
        case 'power_series'
        Expansion=zeros(K*K+2,M,N);
        for n=1:N
            for k1=0:K
                for k2=0:K
                    for m=1:M
                        Expansion(K*k1+k2+1,m,n)=X(NODE,m)^(k1)*X(n,m)^(k2);
                    end
                end
            end
        end
        
    case 'RBF'
        Expansion=zeros(K,M,N);
        for n=1:N
            A=vertcat(X(n,:),X(NODE,:));
            for m1=1:K
                for m2=1:M
                    Expansion(m1,m2,n)=sqrt(2+norm(A(:,m1)-A(:,m2),2)^2);
                end
            end
        end
end

end
