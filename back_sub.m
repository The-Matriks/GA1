function [x] = back_sub(U, b)

  [m, n] = size(U);
  x = zeros(n, 1);
  for i = n:-1:1
    sum=U(i, i+1:n)*x(i+1:n);
    x(i)=(b(i)-sum)/U(i, i);
  end

end
