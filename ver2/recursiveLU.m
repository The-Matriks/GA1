function [L, U] = recursiveLU(A)
    n = size(A, 1);
    if n == 1

        L = 1;
        U = A;
    else

        k = floor(n/2);
        A11 = A(1:k, 1:k);
        A12 = A(1:k, k+1:n);
        A21 = A(k+1:n, 1:k);
        A22 = A(k+1:n, k+1:n);
        

        [L11, U11] = recursiveLU(A11);
        

        L12 = zeros(k, n-k);
        U12 = U11 \ A12;
        

        L21 = A21 / U11;
        U21 = zeros(n-k, k);
        
        [L22, U22] = recursiveLU(A22 - L21 * U12);
        
        L = [L11, L12; L21, L22];
        U = [U11, U12; U21, U22];
    end
end
