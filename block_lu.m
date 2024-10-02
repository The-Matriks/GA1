function [L, U] = block_lu(A, blockSize)
    [n, n] = size(A);
   
    L = eye(n);
    U = zeros(n);
    
    for i = 1:blockSize:n
        iblock = i:min(i+blockSize-1, n);
        
        [L11, U11] = lu(A(iblock, iblock));
        
        L(iblock, iblock) = L11;
        U(iblock, iblock) = U11;
        
        if i + blockSize <= n
            remaining = i+blockSize:n;
            
            L12 = L11 \ A(iblock, remaining);
            U(iblock, remaining) = L12;
            
            U21 = A(remaining, iblock) / U11;
            L(remaining, iblock) = U21;
            
            A(remaining, remaining) = A(remaining, remaining) - U21 * L12;
        end
    end
end
