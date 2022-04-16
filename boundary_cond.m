function [n_, u_] = boundary_cond(n,u,d,dx,dy,a,t,dt)
    
%     n_ = n;
%     u_ = u;
    
    [n_, u_] = reflective_boundary(n,u,d,dx,dy);
%     [n_, u_] = wavemaker_boundary(n,u,d,dx,dy,a,t,dt);
%     n_ = constant_boundary(n,k);

end

