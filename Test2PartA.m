clear all; close all; clc;
%% Problem A2

%% Initial data etc.
N = 12;
a = -1; % left end point of interval
b = 1; % right
TOL = 1e-3;
delta = 0.1;

%% The main loop
while N < 1e+4
    h = abs(b-a)/N;
    x = a:h:b; % node coords
    A=Stiffness_assembler(x);
    B=Load_vector(x);
    M=Mass_assembler(x);
    xi = A\B; % solve system of equation
    zeta = -inv(M)*A*xi; % the discrete laplacian
    eta = zeros(N ,1); % allocate element residuals
    eta2 = zeros(N,1); % allocate element residuals
    for i = 1:N % loop over elements
        h = x(i+1) - x(i); % element length
        if 10*x(i)*sin(7*pi*x(i)) > abs(x(i))
            f = @(x) abs(x)/delta; 
        elseif 10*x*sin(7*pi*x(i)) < abs(x(i))
            f = @(x) -abs(x)/delta;
        else
            f = @(x) (10*x*sin(7*pi*x))/delta;
        end
        R(i) = f(x(i)) + delta*zeta(i);
        a2 = f(x(i)) + delta*zeta(i); % temporary variables
        b2 = f(x(i+1)) + delta*zeta(i+1);
        t = integral(@(x) a2^2 + b2^2, x(i), x(i+1)); % integrate f^2. Trapezoidal rule
        eta(i) = h*sqrt(t);
        eta2(i) = h^2*t; % element residual
    end
    sum_eta2(N) = sum(eta2, 'all');
    n_nodes(N) = length(x);
    if sum_eta2(N) < TOL
        break
    end
    alpha = 0.9; % refinement parameter
    for i = length(eta2)
        if eta2(i) > alpha*max(eta2) % if large residual
            x = [x (x(i+1)+x(i))/2]; % insert new node point
        end 
    end
    x = sort(x); % sort node points accendingly
    N = N+1;
end

%% PLOTS
% First figure
figure(1)
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
grid on
axis padded

subplot(2,2,1)
plot(x,xi)
title('The solution $u_{h}$', 'Interpreter', 'latex')
xlabel('x','Interpreter','latex');ylabel('$u_{h}$','Interpreter','latex');

subplot(2,2,2)
plot(x(2:end), R, 'r')
title('The residual $R(u_{h}) = f + \delta \Delta_{h} u_{h}$', 'Interpreter', 'latex')
xlabel('x','Interpreter','latex');ylabel('$R(u_{h})$','Interpreter','latex');

subplot(2,2,3)
plot(x(2:end),eta, 'g')
title('The error indicator $ \eta (u_{h}) $', 'Interpreter','latex')
xlabel('x','Interpreter','latex');ylabel('$\eta (u_{h})$','Interpreter','latex');

subplot(2,2,4)
plot(x(2:end),round([1./diff(x)]), 'black')
title('The grid size distribution', 'Interpreter','latex')
xlabel('x','Interpreter','latex');ylabel('Grid size distribution','Interpreter','latex');

% Next figure
figure(2)
loglog(n_nodes, sum_eta2 , 'r')
hold on
loglog(n_nodes, n_nodes.^(-1), 'b')
title('The convergence rate', 'Interpreter','latex')
set(gcf,'color','w');
legend('$\sum_{i=1}^{n} \eta_{i}^{2}$', '$\frac{1}{n}$','Interpreter','latex','Location','best');
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of nodes n','Interpreter','latex');
grid on
axis padded