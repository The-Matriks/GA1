function [Z, V_k, eigenvalues, k, cumulative_explained_variance, X_std] = pca(X, threshold)
  [X_std, sigma] = standardize_and_cov(X);

  [eigenvalues, eigenvectors] = qr_algorithm(sigma, 1e-6, 100);

  [eigenvalues_sorted, idx] = sort(eigenvalues, 'descend'); 
  eigenvectors_sorted = eigenvectors(:, idx); 

  explained_variance = eigenvalues_sorted / sum(eigenvalues_sorted);
  cumulative_explained_variance = cumsum(explained_variance);

  k = find(cumulative_explained_variance >= threshold, 1, 'first');

  if isempty(k)
    k = length(eigenvalues_sorted);
  end

  V_k = eigenvectors_sorted(:, 1:k);
  Z = X_std * V_k;
end
