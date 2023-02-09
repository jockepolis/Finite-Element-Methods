clc; clear all; close all;      
%% Problem B2

%% Initial data etc.
fig1 = figure;
geometry = @circleg;
hmax = {1/5, 1/20, 1/40};
T = 2;
delta_1 = 0.01;
alpha = 4;
k = 0.005;
%% Functions
S = @(xi) xi.^2 + xi./(xi + alpha);
a = @(x,y) delta_1;
w = @(x,y) 1; % constant predators
%% Mesh, matrices and solution
for n=1:numel(hmax)
    [p,e,t] = initmesh(geometry, 'hmax', hmax{n});
    A = Stiffness_Assembler_2D(p, t, a);
    M = Mass_assembler_2D(p,t);
    nt = size(t,2);
    np = size(p,2);
    w = rand(np, 1);
    v = 20.*w;
    xi = 1+v; % initial value
    f = @(x,y) xi;
    for l=1:round(T/k)
        for K = 1:nt
            loc2glb = t(1:3,K);
            x = p(1,loc2glb);
            y = p(2,loc2glb); 
            area = polyarea(x,y); 
            bK(K) = 2*sum([f(x(1),y(1));
                           f(x(2),y(2));
                           f(x(3),y(3))])*(area/3);
        end
        integral(l) = mean(bK);
        LHS = 2*M + k*A - k*M;
        RHS = 2*M*xi - k*A*xi + k*M*xi - 2*k*M*S(xi);
        sol = LHS\RHS;
        f = @(x,y) sol;
        xi = sol(1:np);
        pdesurf(p,t,xi) % plot
        pause(0.1)
    end
    if n == 1
        integral1 = integral;
    elseif n == 2
        integral_2 = integral;
    else
        integral3 = integral;
    end
end
title('The solution $v_h$','Interpreter','latex');
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');
%% Integral calculations and 2nd figure
fig2 = figure;
time = linspace(0,2,400);
plot(time, integral1, 'r', LineWidth=2) %integral + pi for predators to be included
hold on
plot(time, integral_2, 'g', LineWidth=2)
hold on
plot(time, integral3, 'b', LineWidth=2)
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
title('The solution $u_{h}$', 'Interpreter', 'latex')
xlabel('x','Interpreter','latex');ylabel('$u_{h}$','Interpreter','latex');
legend('Population rate for $h_{max} = \frac{1}{5}$', 'Population rate for $h_{max} = \frac{1}{20}$', ...
    'Population rate for $h_{max} = \frac{1}{40}$', 'Interpreter','latex','Location','best');
grid on
axis padded