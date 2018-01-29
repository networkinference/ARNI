function dydt = roessler (t,y)

dydt = zeros(size(y));

J=dlmread('Data/connectivity.dat');
[N,~]=size(J);

for i=1:N
    sum=0.0;
    for j=1:N
        sum=J(i,j)*(sin(y(3*(j-1)+1)))+ sum;
    end
    dydt(3*(i-1)+1)= -y(3*(i-1)+2) -y(3*(i-1)+3) + sum;
    dydt(3*(i-1)+2)= y(3*(i-1)+1)+0.1*y(3*(i-1)+2);
    dydt(3*(i-1)+3)= 0.1 + y(3*(i-1)+3)*(y(3*(i-1)+1)-18);
end
end
