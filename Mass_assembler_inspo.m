function M = Mass_assembler_inspo(x)
    n = length(x) - 1;
    M = sparse(zeros(n+1,n+1));
    for i = 1:n
        h = x(i+1) - x(i);
        M(i,i) = M(i,i) + h/3;
        M(i,i+1) = M(i,i+1) + h/6;
        M(i+1,i) = M(i+1,i) + h/6;
        M(i+1,i+1) = M(i+1,i+1) + h/3;
    end
end