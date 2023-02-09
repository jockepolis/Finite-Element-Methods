clc;clear all;close all;
x = linspace(-1,1,1000);
for i=1:1000
    if 10*x(i)*sin(7*pi*x(i)) > abs(x(i))
        b = @(x) abs(x); 
    elseif 10*x*sin(7*pi*x(i)) < -abs(x(i))
        b = @(x) -abs(x);
    else
        b = @(x) 10*x*sin(7*pi*x);
    end
    f(i) = b(x(i));
end
% plot(x,g)
% hold on
pdepe