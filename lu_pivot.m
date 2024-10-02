function [L, U, p] = lu_pivot(A)

  [m, n]  = size(A);
  p = 1:n;
  L = eye(n);

  for k = 1:n-1
    [maxval, maxindex] = max(abs(A(k:n, k)));
    maxindex = maxindex + k - 1;
    A([k, maxindex], :) = A([maxindex, k], :);
    L([k, maxindex], 1:k-1) = L([maxindex, k], 1:k-1);
    p([k, maxindex]) = p([maxindex, k]);
    for i = k+1:n 
      L(i, k) = A(i, k) / A(k, k);
      A(i, k:n) = A(i, k:n) - L(i, k) * A(k, k:n);
    end
  end
  U = A;
end