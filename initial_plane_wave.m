function n = initial_plane_wave(n,A0)
    dims = size(n);
    
    t = 0:dims(2)-1;
    a = A0;
    b1 = round(dims(2)/3);
    c = 100/dims(2);
    
    for i=1:dims(1)
        n(i,:) = n(i,:) + a * exp(-(t-b1).^2 / (2*c^2));
    end
end

