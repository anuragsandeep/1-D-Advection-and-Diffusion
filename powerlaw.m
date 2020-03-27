%% Power law
%% Anurag Sandeep K. (UIN:624008228)


function [T]=powerlaw(ITMAX,dx,rho,u,gama,tau,d)

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

% power law formulation
% 1. A(|P|)
AP=max(0,(1-0.1*abs(Pe))^5);
APend=max(0,(1-0.1*abs(Peend))^5);

aW=D*AP+max(F,0);
aE=D*AP+max(-F,0);

% Calling Tri-diagonal
A(1,2)=(Dend*APend+max(F,0))+aE + tau*2*dx/d;
A(1,3)=-aE;
B(1)=(Dend*APend+max(F,0))*T(1);

for i=3:ITMAX-2
    A(i-1,1)=-aW;
    A(i-1,2)=aW+aE + tau*2*dx/d;
    A(i-1,3)=-aE;
    B(i-1)=0; 
end

A(ITMAX-2,1)=-aW;
A(ITMAX-2,2)=(Dend*APend+max(-F,0))+aW + tau*2*dx/d;
B(ITMAX-2)=(Dend*APend+max(-F,0))*T(end);

T(2:end-1)=Tridiagonal(ITMAX-2,A,B);

end
