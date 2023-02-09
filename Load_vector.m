function b=Load_vector(x)
delta = 0.1;
N = length(x) - 1; 
b = zeros(N+1, 1);
for i = 1:N
    if 10*x(i)*sin(7*pi*x(i)) > abs(x(i))
        f = @(x) abs(x)/delta; 
    elseif 10*x(i)*sin(7*pi*x(i)) < -abs(x(i))
        f = @(x) -abs(x)/delta;
    else
        f = @(x) (10.*x.*sin(7*pi.*x))/delta;
    end
    h = x(i+1) - x(i);
    n = [i i+1];
    b(n) = b(n) + [f(x(i)); f(x(i+1))]*h/2;
end
end