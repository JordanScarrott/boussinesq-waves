

% d is boundary depth
% a is wave amplitude
% t is current time
function [n_,u_] = wavemaker_boundary(n,u,d,dx,dy,A,t_curr,dt)

    dims_n = size(n);

    n_ = n;
    u_ = u;
%     d = 2;
        
    % Math constants
    h0 = 0.45;
    w = 1;
    k = 1;
    a = -0.390;
    
    c = 1;
    
    
    t = ((d-1:-1:0) * dt + t_curr) .* ones(dims_n(1),d);
    
    % Wave elevation
    n_(:,1:d) = A * sin(c*t);
    
    % The correct wave speed for this elevation
    coeff = (w)/(k*h0*(1-(a+1/3)*(k*h0)^2));
    
    u_.u(:,1:d) = coeff * n_(:,1:d).*cos(c*t);
    u_.v(:,1:d) = coeff * n_(:,1:d).*sin(c*t);
    
            

end

