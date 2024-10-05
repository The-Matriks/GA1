function [Q, R] = householder(A)
  [m, n] = size(A);
  Q = eye(m);
  R = A;

  for j = 1:n
    x = R(j:m, j);
    normx = norm(x);

    if normx == 0
        continue;
    end
    s = -sign(x(1));
    u1 = x(1) - s * normx;
    w = x / u1;
    w(1) = 1;
    w = w / norm(w);

    R(j:m, j:n) = R(j:m, j:n) - 2 * (w * (w' * R(j:m, j:n)));

    Q(:, j:m) = Q(:, j:m) - 2 * (Q(:, j:m) * w) * w';
  end

  Q = Q';
end
