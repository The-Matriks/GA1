function [L, U] = blockLU(A, blockSize)

    n = size(A, 1);
    L = eye(n);
    U = zeros(n);

    for i = 1:blockSize:n
        j = min(i + blockSize - 1, n);
        
        [L(i:j, i:j), U(i:j, i:j)] = lu(A(i:j, i:j), 'vector');
    
        if j < n
            L(j+1:n, i:j) = A(j+1:n, i:j) / U(i:j, i:j);
            U(i:j, j+1:n) = L(i:j, i:j) \ A(i:j, j+1:n);
            A(j+1:n, j+1:n) = A(j+1:n, j+1:n) - L(j+1:n, i:j) * U(i:j, j+1:n);
        end
    end
end
