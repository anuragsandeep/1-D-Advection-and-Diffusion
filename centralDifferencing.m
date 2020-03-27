%% Central differencing scheme
%% Anurag Sandeep K. (UIN:624008228)


function [T]=centralDifferencing(ITMAX,dx,rho,u,gama,tau,d)

T=zeros(ITMAX,1);
T(1)=1000;
T(end)=700;
A=[];B=[];

% Flow strength
F=rho*u;

% Diffusion strength
D=gama/dx;

% Peclet number
Pe=F/D;

% central differencing
aW=D+(F/2);
aE=D-(F/2);

% Calling Tri-diagonal
A(1,2)=2*D+(F/2) + D-(F/2) + tau*2*dx/d;
A(1,3)=-aE;
B(1)=(2*D+(F/2))*T(1);

for i=3:ITMAX-2
    A(i-1,1)=-aW;
    A(i-1,2)=aW+aE + tau*2*dx/d;
    A(i-1,3)=-aE;
    B(i-1)=0; 
end

A(ITMAX-2,1)=-aW;
A(ITMAX-2,2)=D+(F/2) + 2*D-(F/2) + tau*2*dx/d;
B(ITMAX-2)=(2*D-(F/2))*T(end);


T(2:end-1)=Tridiagonal(ITMAX-2,A,B);

end