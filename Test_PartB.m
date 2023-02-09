clc; clear all; close all;
%% Problem B2

%% Initial data etc.
fig1 = figure;
geometry = @circleg;
hmax = {1/5, 1/20, 1/40};
T = 2;
delta_1 = 0.01;
alpha = 4;
k = 0.05;
time = linspace(0,2,round(T/k));
%% Functions
S = @(xi) xi.^2 + xi./(xi + alpha);
a = @(x,y) delta_1;
%% Mesh, matrices and solution
for n=1:numel(hmax)
    [p,e,t] = initmesh(geometry, 'hmax', hmax{n});
    x = p(1,:)'; y = p(2,:)';
%     k = hmax{n}/5;
    A = Stiffness_Assembler_2D(p, t, a);
    M = Mass_assembler_2D(p,t);
    nt = size(t,2);
    np = size(p,2);
    w = rand(np, 1);
    v = 20.*w;
    xi = 1+v; % initial value
    for l=1:round(T/k)
        LHS = 4*M + k*A - k*M;
        RHS = 4*M*xi - k*A*xi + k*M*xi - 4*k*M*S(xi);
        sol = LHS\RHS;
        xi = sol(1:np);
        for K = 1:nt
            loc2glb = t(1:3,K);
            x = p(1,loc2glb);
            y = p(2,loc2glb); 
            area = polyarea(x,y); 
            bK = (sum(xi(loc2glb))*area)/3; 
        end
        integral(l) = bK;
    end
    plot(time, integral)
    hold on
end
title('The solution $v_h$','Interpreter','latex');
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');
%% Integral calculations and 2nd figure
% fig2 = figure;
% % time1 = linspace(0,2,50);
% % plot(time1, integral1, 'r')
% % hold on
% time2 = linspace(0,2,200);
% plot(time2, integral2, 'g')
% hold on
% time3 = linspace(0,2,400);
% plot(time3, integral3, 'b')


