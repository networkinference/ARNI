function dydt = michaelis_menten (t,y)

dydt = zeros(size(y));

J=dlmread('Data/connectivity.dat');

k=length(J);

for i=1:k
    sum=0.0;
    for j=1:k
        sum=J(i,j)*(y(j)/(1+y(j))) + sum;
        %sum=J(i,j)*sin(y(j)-y(i)) + sum;
        
    end
    dydt(i)=-y(i)+ sum;
    %dydt(i)=w(i)+ tanh(sum);
    
end
end
