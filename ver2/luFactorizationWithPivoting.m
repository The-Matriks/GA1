function [L, U, P] = luFactorizationWithPivoting(A)
    n = size(A, 1);  
    L = eye(n);      
    U = A;           
    P = eye(n);      
    for k = 1:n-1

        [~, m] = max(abs(U(k:n, k)));
        m = m + k - 1;
        if m ~= k

            U([k m], :) = U([m k], :);
            P([k m], :) = P([m k], :);
            
            if k >= 2
                L([k m], 1:k-1) = L([m k], 1:k-1);
            end
        end
        

        for j = k+1:n
            L(j, k) = U(j, k) / U(k, k);
            U(j, k:n) = U(j, k:n) - L(j, k)*U(k, k:n);
        end
    end
    

    L(:, n) = 0;
    L(n, n) = 1;
end
