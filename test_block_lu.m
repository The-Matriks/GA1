A1 = [4, 3; 6, 3];
A2 = [1, 2, 3; 4, 5, 6; 7, 8, 10];
A3 = [1, 2, 3; 4, 5, 6; 7, 8, 9];

[L1, U1] = block_lu(A1, 2);
[L2, U2] = block_lu(A2, 3);

disp('Verification for A1:');
disp('Original A1:');
disp(A1);
disp('Computed LU for A1:');
disp(L1 * U1);

disp('Verification for A2:');
disp('Original A2:');
disp(A2);
disp('Computed LU for A2:');
disp(L2 * U2);

try
    [L3, U3] = block_lu(A3, 3);
    disp('Verification for A3:');
    disp('Original A3:');
    disp(A3);
    disp('Computed LU for A3:');
    disp(L3 * U3);
catch exception
    disp('Error encountered for A3:');
    disp(exception.message);
end