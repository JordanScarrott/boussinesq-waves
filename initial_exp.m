function n = initial_exp(n, A0, x, y)
    dims = size(n);
    
    X = x .* ones(dims(1), dims(2));
    Y = y' .* ones(dims(1), dims(2));
    
    n = A0 * exp(-2*((X - 3.75).^2 + (Y - 3.75).^2));
    
end

