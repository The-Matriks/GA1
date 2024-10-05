warning('off', 'all');

% N = [16];
N = [16 64 128 256 512 1000];

for i=1:length(N)
    n = N(i);
    P = [1 2 3 n/4 n/2];
    Q = [1 1 4 n/4 n/2];
    
    for j=1:length(P)
        printf('n=%d, p=%d, q=%d\n', n, P(j), Q(j));
        p = P(j);
        q = Q(j);

        matrix = create_band_matrix(n, p, q);
        while det(matrix) == 0 || isinf(cond(matrix))
            printf('det(matrix) = 0, cond(matrix) = %f\n', cond(matrix));
            matrix = create_band_matrix(n, p, q);
        end
        b = rand(n, 1);

        % LU Pivot + Forward
        start = tic();
        [L, U, p] = lu_pivot(matrix);

        x = forward_elim(L, b);
        lu_pivot_forward_time = toc(start);

        % LU Pivot + Backward
        start = tic();
        [L, U, p] = lu_pivot(matrix);

        y = back_sub(U, b);
        lu_pivot_backward_time = toc(start);

        disp('------------------------------------------------------------');
        printf('LU Pivot Error ');
        disp(norm(p*matrix-L*U));
        printf('LU Pivot + Forward Error ');
        disp(norm(L*x-b));
        printf('LU Pivot + Backward Error ');
        disp(norm(U*y-b));
        printf('LU Pivot + Forward Time ');
        disp(lu_pivot_forward_time);
        printf('LU Pivot + Backward Time ');
        disp(lu_pivot_backward_time);
        disp('------------------------------------------------------------');

        % Recursive LU + Forward
        start = tic();
        [L, U] = recursive_lu(matrix);

        x = forward_elim(L, b);
        recursive_lu_forward_time = toc(start);

        % Recursive LU + Backward
        start = tic();
        [L, U] = recursive_lu(matrix);

        y = back_sub(U, b);
        recursive_lu_backward_time = toc(start);

        disp('------------------------------------------------------------');
        printf('Recursive LU + Forward Error ');
        disp(norm(L*x-b));
        printf('Recursive LU + Backward Error ');
        disp(norm(U*y-b));
        printf('Recursive LU + Forward Time ');
        disp(recursive_lu_forward_time);
        printf('Recursive LU + Backward Time ');
        disp(recursive_lu_backward_time);
        disp('------------------------------------------------------------');

        % Block LU + Forward
        start = tic();
        [L, U] = block_lu(matrix, size(matrix, 1), (n/2)*(n/2));

        x = forward_elim(L, b);
        block_lu_forward_time = toc(start);

        % Block LU + Backward
        start = tic();
        [L, U] = block_lu(matrix, size(matrix, 1), (n/2)*(n/2));

        y = back_sub(U, b);
        block_lu_backward_time = toc(start);

        disp('------------------------------------------------------------');
        printf('Block LU + Forward Error ');
        disp(norm(L*x-b));
        printf('Block LU + Backward Error ');
        disp(norm(U*y-b));
        printf('Block LU + Forward Time ');
        disp(block_lu_forward_time);
        printf('Block LU + Backward Time ');
        disp(block_lu_backward_time);
        disp('------------------------------------------------------------');
    end
end
