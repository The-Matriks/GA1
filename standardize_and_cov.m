function [X_std, sigma] = standardize_and_cov(X)
    % This function standardizes each column of X and computes the covariance matrix.
    % Input: 
    %   X - m x n data matrix (m samples, n features)
    % Output: 
    %   X_std - standardized version of X (mean 0, variance 1 for each column)
    %   sigma - covariance matrix of the standardized data

    [m, n] = size(X);
    
    % Step 1: Compute the mean and standard deviation of each column
    mu = mean(X, 1);            % 1xN row vector of means for each column
    sigma_j = std(X, 0, 1);      % 1xN row vector of std devs (ddof=0 for population std)

    % Step 2: Standardize the data (for each column: subtract mean and divide by std dev)
    X_std = (X - mu) ./ sigma_j; % Standardize each column

    % Step 3: Compute the covariance matrix of the standardized data
    sigma = (X_std' * X_std) / (m - 1);  % Covariance matrix
end
