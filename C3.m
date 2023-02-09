
clear all; close all; clc;
a = 0; b = 2*pi;
T = 2; %Final time
BETA = @(u) (1/2.* u .^2);

colors=[28/255 104/255 192/255;
    184/255 19/255 32/255;
    71/255 71/255 71/255;
    234/255 97/255 25/255;
    255/255 99/255 248/255;
    66/255 112/255 28/255;
    2/255 3/255 9/255];

g=0;
fig1=figure;
hold on
for n = [41, 81, 161, 321, 641]
    g=g+1;
    h = abs(a-b)/(n-1); % Spatialsteps
    x = linspace(a,b,n); % Space-vector
    eps = 0.5*h;
    C = 0.1;
    k = C * h^2; %Timestep
    t=0:k:T; %Time-vector
    
    M = Mass_assembler(x); %Matrix assembly
    A = Advection_assembler(x);
    S = StiffnessAssembler1D(x);
    
    U=zeros(length(x),length(t)); %Preallocating the solution matrix
    U(:,1) = sin(x); %inserting initial data
    
    for i = 1:length(t)
        %RK4 with the Conjugate-Gradient method ("cg")
        K1 = cg(M,(-A*BETA(U(:,i)) - eps*S*U(:,i)),1e-3);
        K2 = cg(M,(-A*BETA(U(:,i)+k/2*K1) - eps*S*(U(:,i)+k/2*K1)),1e-3);
        K3 = cg(M,(-A*BETA(U(:,i)+k/2*K2) - eps*S*(U(:,i)+k/2*K2)),1e-3);
        K4 = cg(M,(-A* BETA((U(:,i)+k*K3)) - eps*S*(U(:,i)+k*K3)),1e-3);
        U(:,i+1)=U(:,i) + k/6*(K1+2*K2+2*K3+K4);
        U(:,i+1)=U(:,i) + k/6*(K1+2*K2+2*K3+K4);
        
        if i<length(t) %Adding the boundary condition for every timestep
            U(1,i+1) = 0;
            U(end,i+1) = 0;
        else
            
        end
    end
    plot(x,U(:,end),'color',colors(g,:))
end

set(gcf,'color','w');
legend('N=41','N=81','N=161',...
    'N=321','N=641',...
    'Interpreter','latex','Location','best');
xlabel('x','Interpreter','latex');ylabel('$u_h(x,T_{final})$','Interpreter','latex');
title(['Solution at final time $T = $ ',num2str(T)],'Interpreter','latex');
xlim([a b]); xticks(a:pi/2:b); xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
set(gca,'TickLabelInterpreter','latex')
grid on
axis padded
hold off

