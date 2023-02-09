clc; clear all; close all;

f = @(x,y) 8*pi^2*sin(2*pi*x)*sin(2*pi*y);
g = @(x,y) sin(2*pi*x)*sin(2*pi*y);
a = @(x,y) 1;

geometry = @circleg;
hmax = 1/2;
[p,e,t] = initmesh(geometry, 'hmax', hmax);
I = eye(length(p));
A = Stiffness_Assembler_2D(p, t, a);
b = Load_vector_2D(p, t, f);
np = size(p,2); % total number of nodes
fixed = unique([e(1,:) e(2,:)]); % boundary nodes
free = setdiff([1:np],fixed); % interior nodes
b = b(free)-A(free,fixed)*g; % modify b
A = A(free,free); % modify A
xi = zeros(np,1); % allocate solution 
xi(fixed) = g; % insert fixed node values 
xi(free) = A\b;
