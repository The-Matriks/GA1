warning('off', 'all');

% Define matrix sizes
N = [16, 64, 128, 256, 512, 1000];

% Loop through each matrix size
for i = 1:length(N)
    n = N(i);
    P = [1, 2, 3, n/4, n/2];
    Q = [1, 1, 4, n/4, n/2];
    
    for j = 1:length(P)
        p = P(j);
        q = Q(j);
        fprintf('n=%d, p=%d, q=%d\n', n, p, q);

        % Generate a matrix and ensure it is non-singular
        matrix = create_band_matrix(n, p, q);
        while det(matrix) == 0 || isinf(cond(matrix))
            fprintf('det(matrix) = 0, cond(matrix) = %f\n', cond(matrix));
            matrix = create_band_matrix(n, p, q);
        end
        
        cond_number = cond(matrix);
        fprintf('Condition Number: %f\n', cond_number);
        b = rand(n, 1);  % Random b vector

        % LU Pivot Decomposition
        [L, U, perm, flops_lu] = lu_pivot(matrix);
        [x, flops_forward] = forward_elim(L, b);
        [y, flops_backward] = back_sub(U, x);
        
        fprintf('LU Pivot Error %e\n', norm(perm * matrix - L * U));
        fprintf('LU Pivot + Forward Error %e\n', norm(L * x - b));
        fprintf('LU Pivot + Backward Error %e\n', norm(U * y - b));
        fprintf('LU Pivot FLOPs Count: %d\n', flops_lu + flops_forward + flops_backward);
        fprintf('LU Pivot FLOPs Complexity: %e\n', 2*n^3/3 + 2*n*(n - 1));
        fprintf('------------------------------------------------------------\n');

        % Recursive LU Decomposition
        try
            [L, U, flops_recursive] = recursive_lu(matrix);
            [x, flops_forward] = forward_elim(L, b);
            [y, flops_backward] = back_sub(U, x);
            
            fprintf('Recursive LU + Forward Error %e\n', norm(L * x - b));
            fprintf('Recursive LU + Backward Error %e\n', norm(U * y - b));
            fprintf('Recursive LU FLOPs Count: %d\n', flops_recursive + flops_forward + flops_backward);
            fprintf('Recursive LU FLOPs Complexity: %e\n', 2*n^3/3 + 2*n*(n - 1));
            fprintf('------------------------------------------------------------\n');
        catch
            fprintf('Recursive LU failed due to matrix issues.\n');
        end

        % Block LU Decomposition
        try
            [L, U, flops_block] = block_lu(matrix, n, max(1, floor(n/4)));
            [x, flops_forward] = forward_elim(L, b);
            [y, flops_backward] = back_sub(U, x);
            
            fprintf('Block LU + Forward Error %e\n', norm(L * x - b));
            fprintf('Block LU + Backward Error %e\n', norm(U * y - b));
            fprintf('Block LU FLOPs Count: %d\n', flops_block + flops_forward + flops_backward);
            fprintf('Block LU FLOPs Complexity: %e\n', 2*n^3/3 + 2*n*(n - 1));
            fprintf('------------------------------------------------------------\n');
        catch
            fprintf('Block LU failed due to matrix issues.\n');
        end
    end
end
