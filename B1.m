clc; clear all; close all;
%% Problem B1

%% Functions
f = @(x,y) 8*pi^2.*sin(2*pi.*x).*sin(2*pi.*y);
g = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
a = @(x,y) 1;
%% Mesh, matrices and solution
geometry = @circleg;
hmax = {1/2, 1/4, 1/8, 1/16, 1/32};
for k=1:numel(hmax)
    [p,e,t] = initmesh(geometry, 'hmax', hmax{k});
    I = eye(length(p));
    A = Stiffness_Assembler_2D(p, t, a);
    b = Load_vector_2D(p, t, f);
    x = p(1,:);
    y = p(2,:);
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:)) = g(p(1,e(1,:)),p(2,e(1,:)));
    xi = A\b;
    u_exact = g(x,y);
    err = u_exact' - xi;
    EnE(k)=sqrt(err'*A*err);
end
%% Plots and calculation of convergence rate
fig1=figure;
hmax_plot = [1/2, 1/4, 1/8, 1/16, 1/32];
p_convergence=log(EnE(end)/EnE(end-1))/log(hmax_plot(end)/hmax_plot(end-1));
hpmax_plot = [1/2, 1/4, 1/8, 1/16, 1/32].^p_convergence;
loglog(hmax_plot,EnE,'ok-','MarkerFaceColor',[0.635 0.078 0.1840]);
hold on
loglog(hmax_plot, hpmax_plot, 'ok-', 'MarkerFaceColor', 'k');
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
xlabel('$h_{max}$','Interpreter','latex');
legend('$|| u_{exact}-u_h ||_2$', '$h_{max}^{p}$','Interpreter','latex','Location','best');
title('$|| u_{exact}-u_h ||_2$ and $h_{max}^{p}$ plotted against the spatial step h','Interpreter','latex');
grid on
axis padded

fig2=figure;
pdeplot(p,e,t,"XYData",xi)
set(gcf,'color','w');
title('The solution $u_h$','Interpreter','latex');
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');

fig3=figure;
pdeplot(p,e,t,"XYData",u_exact)
set(gcf,'color','w');
title('The solution $u_{exact}$','Interpreter','latex');
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');