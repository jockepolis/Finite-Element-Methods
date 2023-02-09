function A = Advection_assembler_Inspo(x)
   n = length(x) - 1;
   A = sparse(zeros(n+1, n+1));
   for i = 1:n
       N = [i i+1];
       A(N,N) = A(N,N) + [0 1/2; -1/2 0];
   end
  A(1,1) = -1/2; A(end,end) = 1/2;
end