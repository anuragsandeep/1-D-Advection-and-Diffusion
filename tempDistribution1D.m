%% Temperature distribution in 1-D with advection and diffusion 
%% with TDMA and two-node formulation
%% Anurag Sandeep K. (UIN:624008228)

clear all 
clc

% INPUT PARAMETERS
L=2;       % length of channel
d=0.05;    % height of channel
rho=10000; % density
cp=140;    % specific heat 
k=21;      % heat transfer co-efficient [W/m/K]
h=2000;    % convective heat transfer co-efficient
u=0.2; % constant velocity
gama=k/cp; % diffusion co-efficient
tau=h/cp;
nCV=9;

% flags for schemes
% 1: central differencing
% 2: exponential
% 3: powerlaw
flag=3;

if (flag==1)
        fprintf('Computing with central differencing scheme\n\n');
    
    elseif (flag==2)
        fprintf('Computing with exponential scheme\n\n');

    else
        fprintf('Computing with power-law\n\n');
    end

count=0;
diff=1;
while (diff > 1e-3)
    count=count+1;
    ITMAX(count)=nCV+2;
    dx(count)=L/(ITMAX(count)-2); % grid-size in x
    
    if (flag==1)
        T{count}=centralDifferencing(ITMAX(count),dx(count),rho,u,gama,tau,d);
    
    elseif (flag==2)
        T{count}=exponential(ITMAX(count),dx(count),rho,u,gama,tau,d);

    else
        T{count}=powerlaw(ITMAX(count),dx(count),rho,u,gama,tau,d);
    end
    
    % Picking the mid-cell to compute the temperature difference
    % between two grids
    if (count>1)

        Tnew=T{count}((ITMAX(count)+1)/2);
        Told=T{count-1}((ITMAX(count-1)+1)/2);
  
    % computing the max temperature difference
    diff=abs(Tnew-Told);
    
    end  
    
    % doubling the number of CV's
    nCV=nCV*2 + 1;
end % while loop

% PLOTTING RESULTS
% constructing the grid locations
for i=1:length(ITMAX)-1
n{i}=zeros(ITMAX(i),1);
n{i}(2)=dx(i)/2;
n{i}(end)=2;
    for j=3:length(n{i})-2
        n{i}(j)=n{i}(j-1)+dx(i);
    end
n{i}(end-1)=n{i}(end)-dx(i)/2;
end

% plotting temperature based on power-law
if (flag==3)
    for i=1:length(ITMAX)-1
    txt = ['CVs = ',num2str(ITMAX(i))];    
    plot(n{i},T{i},'MarkerSize',3,'Marker','o','DisplayName',txt); hold on
    end  
    xlabel('location (m)')
    ylabel('Temperature (^{\circ}C)','Interpreter','tex')
    title('Convergence for Temperature based on Power-law')
end
legend show












