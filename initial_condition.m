

function n = initial_condition(n,A0,x,y)
%     n = initial_gaussian_2d(n,A0);
%     n = initial_plane_wave(n,A0);
%     n = initial_sech(n,A0,0);
    n = initial_exp(n,A0,x,y);
end

