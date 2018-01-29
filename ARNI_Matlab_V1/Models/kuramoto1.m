function dydt = kuramoto1 (t,y)

dydt = zeros(size(y));

w=dlmread('Data/frequencies.dat');

k=length(w);

J=dlmread('Data/connectivity.dat');


for i=1:k
    sum=0.0;
    for j=1:k
        sum=J(i,j)*sin(y(j)-y(i)) + sum;
    end
    dydt(i)=w(i)+ sum;
    
end
end
