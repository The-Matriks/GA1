function [L, U] = block_lu(A, n, r)
    if n <= r
        L = eye(n);
        U = zeros(n);
        for i = 1:n
            for j = i:n  
                sum = 0;
                for k = 1:i-1
                    sum = sum + L(i,k) * U(k,j);
                end
                U(i,j) = A(i,j) - sum;
            end
            for j = i+1:n 
                sum = 0;
                for k = 1:i-1
                    sum = sum + L(j,k) * U(k,i);
                end
                L(j,i) = (A(j,i) - sum) / U(i,i);
            end
        end
    else

        [L11, U11] = block_lu(A(1:r, 1:r), r, r);
        

        U12 = L11 \ A(1:r, r+1:end);
        L21 = A(r+1:end, 1:r) / U11;

        A22 = A(r+1:end, r+1:end) - L21 * U12;
        

        [L22, U22] = block_lu(A22, n-r, r);

        L = [L11, zeros(r, n-r); L21, L22];
        U = [U11, U12; zeros(n-r, r), U22];
    end
end
