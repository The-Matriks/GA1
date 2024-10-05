import numpy as np
import pandas as pd
from numpy.linalg import norm

def standardize(X):
    if X.size == 0:
        raise ValueError("Input array X is empty.")
    mean = np.mean(X, axis=0)
    std_dev = np.std(X, axis=0)
    if np.any(std_dev == 0):
        raise ValueError("Standard deviation is zero for one or more features, cannot standardize.")
    X_std = (X - mean) / std_dev
    return X_std

def covariance_matrix(X):
    if X.size == 0:
        raise ValueError("Input array X is empty.")
    n_samples = X.shape[0]
    return (1 / (n_samples - 1)) * X.T @ X

def qr_decomposition_householder(A, max_iter=100, tol=1e-8):
    n = A.shape[0]
    Q = np.eye(n)
    R = A.copy()
    for i in range(min(n, max_iter)):
        x = R[i:, i]
        if np.max(np.abs(x)) < tol:
            break
        e = np.zeros_like(x)
        e[0] = norm(x)
        u = x - e
        u /= norm(u)
        H = np.eye(n)
        H[i:, i:] -= 2.0 * np.outer(u, u)
        R = H @ R
        Q = Q @ H.T
    return Q, R

def perform_pca(X, k=6):
    X_std = standardize(X)
    cov_matrix = covariance_matrix(X_std)
    eigenvalues, eigenvectors = qr_decomposition_householder(cov_matrix)
    eigenvalues = np.diag(eigenvalues)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    V_k = eigenvectors[:, :k]
    Z = X_std @ V_k
    return Z, k, eigenvalues, eigenvectors

def read_csv_data(file_path):
    df = pd.read_csv(file_path)
    X = df.select_dtypes(include=[np.number]).values
    return X

csv_file = "pokindex_data.csv"
df = read_csv_data(csv_file)
X = df[:, :-1]

Z, k, eigenvalues, eigenvectors = perform_pca(X)

print("Number of selected components (k):", k)
print("Transformed data shape:", Z.shape)

df = pd.read_csv("pokindex_data.csv")