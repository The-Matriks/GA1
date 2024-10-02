function [ B ] = create_band_matrix(n, p, q)
    B = zeros(n);
    for i=-q:p
        B = B + diag(rand(n-abs(i),1),-i);
    end
end
