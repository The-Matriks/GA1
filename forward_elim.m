function [x] = forward_elim(L, b)

  [m, n] = size(L);
  x = zeros(n, 1);

  for i = 1:n
    sum=L(i, 1:i-1)*x(1:i-1);
    x(i)=(b(i)-sum)/L(i, i);
  end

end