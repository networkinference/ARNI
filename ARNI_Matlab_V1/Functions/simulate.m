function simulate(MODEL,N,NI,S,M)
% simulate(MODEL,N,NI,S,M) generates time series of networks of dynamical
% systems for several different intial conditions.
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
%
% Input type
% ------------------
% MODEL: string
% N:     integer
% NI:    integer (NI<N)
% S:     integer
% M:     integer
%
% Output
% ------------------
% 'Data/data.dat':     File containing all simulated time series in a
%                      concatenaded form.
% 'Data/ts_param.dat': File containing time series parameters, i.e. S and
%                      M, for later extracting the different time series.
%
% Example
% ------------------
% simulate('kuramoto2',25,4,30,10) generates 30 time series of 10 time
% points each for a network of 25 oscillators defined by the model
% kuramoto2. Each oscillator has 4 incoming connections.
%
% Accompanying material to "Model-free inference of direct interactions 
% from nonlinear collective dynamics".
%
% Author: Jose Casadiego
% Date:   May 2017

% Sampling rate of time series
resolution=1;

models={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
if any(ismember(models,MODEL))==0
    disp('ERROR: MODEL must be a valid string: kuramoto1, kuramoto2, michaelis_menten, roessler');
else
    system('rm -r Data');
    system('mkdir Data');
    disp('Creating network structure...')
    topology(N,'homogeneous','directed',NI);
    disp('Simulating time series...')
    Y=[];
    switch MODEL
        case 'kuramoto1'
            w=-2 + (4).*rand(N,1);
            dlmwrite('Data/frequencies.dat', w, 'delimiter', '\t', 'precision', 4);
            for s=1:S
                init=-3.14 +(3.14+3.14)*rand(N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('kuramoto1',tspan,init);
                Y=[Y;y];
            end
        case 'kuramoto2'
            w=-2 + (4).*rand(N,1);
            dlmwrite('Data/frequencies.dat', w, 'delimiter', '\t', 'precision', 4);
            for s=1:S
                init=-3.14 +(3.14+3.14)*rand(N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('kuramoto2',tspan,init);
                Y=[Y;y];
            end
        case 'michaelis_menten'
            for s=1:S
                init=1+1*rand(N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('michaelis_menten',tspan,init);
                Y=[Y;y];
            end
        case 'roessler'
            for s=1:S
                init=-5 +(5+5)*rand(3*N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('roessler',tspan,init);
                Y=[Y;y];
            end     
    end
    ts_param=[S,M];
    dlmwrite('Data/data.dat', Y, 'delimiter', '\t');
    dlmwrite('Data/ts_param.dat', ts_param, 'delimiter', '\t');
    clear;
    disp('Simulation finished!');
end
end