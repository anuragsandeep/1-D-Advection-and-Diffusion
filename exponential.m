%% Exponential scheme/exact solution
%% Anurag Sandeep K. (UIN:624008228)


function [T]=exponential(ITMAX,dx,rho,u,gama,tau,d)

T=zeros(ITMAX,1);
T(1)=1000;
T(end)=700;
A=[];B=[];

% Flow strength
F=rho*u;

% Diffusion strength
D=gama/dx;
Dend=2*gama/dx;

% Peclet number
Pe=F/D;
Peend=F/Dend;

% exponential scheme
aW=F+F/(exp(Pe)-1);
aE=F/(exp(Pe)-1);
aP=aW+aE;

% Calling Tri-diagonal
A(1,2)=(F+F/(exp(Peend)-1)) + F/(exp(Pe)-1) + tau*2*dx/d;
A(1,3)=-aE;
B(1)=(F+F/(exp(Peend)-1))*T(1);

for i=3:ITMAX-2
    A(i-1,1)=-aW;
    A(i-1,2)=aW+aE + tau*2*dx/d;
    A(i-1,3)=-aE;
    B(i-1)=0; 
end

A(ITMAX-2,1)=-aW;
A(ITMAX-2,2)=(F+F/(exp(Pe)-1)) + F/(exp(Peend)-1) + tau*2*dx/d;
B(ITMAX-2)=F/(exp(Peend)-1)*T(end);


T(2:end-1)=Tridiagonal(ITMAX-2,A,B);

end