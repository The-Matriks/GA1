function A = generateBandedMatrix(N, p, q)
    A = zeros(N);  

    A = A + diag(10*rand(N,1) + 10); 

    for i = 1:p  
        A = A + diag(0.5*rand(N-i,1), -i);
    end
    for i = 1:q  
        A = A + diag(0.5*rand(N-i,1), i);
    end
end
