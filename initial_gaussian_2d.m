function n = initial_gaussian_2d(n, A0)
    dims = size(n);
    a = A0;
    b1 = round(dims(2)/2);
    b2 = round(dims(1)/2);
    c = sqrt(dims(1));
    
    for k=1:dims(1)
        for j=1:dims(2)
            n(k,j) = n(k,j) + a * exp(-((j-b1).^2 + (k-b2).^2) / (2*c^2));
        end
    end
end

