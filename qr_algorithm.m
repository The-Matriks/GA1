function [eigenvalues, eigenvectors] = qr_algorithm(A, tol, max_iter)
    [m, n] = size(A);
    Q_total = eye(m);
    Ak = A; 
    
    for k = 1:max_iter
        [Q, R] = householder(Ak);
        Ak = R * Q; 
        Q_total = Q_total * Q;
        
        if norm(Ak - diag(diag(Ak)), 'fro') < tol
            break;
        end
    end
    
    eigenvalues = diag(Ak);
    eigenvectors = Q_total;
    for i = 1:size(eigenvectors, 2)
        eigenvectors(:, i) = eigenvectors(:, i) / norm(eigenvectors(:, i));
    end
end
