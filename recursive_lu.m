function [L, U] = recursive_lu(A)
    [n, n] = size(A);
    
    if n <= 2
        [L, U] = lu(A);
        return;
    end

    mid = floor(n/2);
    A11 = A(1:mid, 1:mid);
    A12 = A(1:mid, mid+1:end);
    A21 = A(mid+1:end, 1:mid);
    A22 = A(mid+1:end, mid+1:end);
    
    [L11, U11] = recursiveLU(A11);

    A12 = L11 \ A12; 
    A21 = A21 / U11; 
    
    A22 = A22 - A21 * A12;

    [L22, U22] = recursiveLU(A22);

    L = [L11, zeros(mid, n-mid); A21, L22];
    U = [U11, A12; zeros(n-mid, mid), U22];
end
