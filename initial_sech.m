
% Remember to set the constants in this function according to Wei1995 when
% doing actual final tests.
% Set is_velocity to 1 if you want to get the velocity function instead
function n = initial_sech(n, A0, is_velocity)
    dims = size(n);
    
    x = 0:dims(2)-1;
    a1 = A0/2;
    a2 = A0/2;
    B = 0.2;
    b1 = dims(2)/4;
    c = 1;
    
    smaller_dim = min(dims(2),dims(1));
    
    if is_velocity == 0
        for i=1:smaller_dim
            n(i,:) = n(i,:) + a1 * sech(B*(c*x-b1)).^2 + a2 * sech(B*(c*x-b1)).^4;
        end
    elseif(is_velocity == 1)
        for i=1:smaller_dim
            n(i,:) = n(i,:) + a1 * sech(B*(c*x-b1)).^2;
        end
    end
end

