

% Handles the boundary condition for reflective boundaries
% They are computed per side with the following equations
% u.n = 0
% grad n.normal_vector = 0
% duT/dn = 0
%
% d is the depth of the boundary (how many cells deep is the boundary)



% Not sure of the first and last boundary conditions. For the first BC I
% swapped all the u.us with u.vs and vice versa
function [n_,u_] = reflective_boundary(n,u,d,dx,dy)
    
    n_dims = size(n);

    n_ = n;
    u_ = u;

    d = 1;
    b = [48; -36; 16; -3];
    p = length(b);

    % Computed per side (left,right,up,down)
    % Left side
    n_(d:end-d+1,1) = (n(d:end-d+1,2:p+1) * b)/25;
    u_.v(d:end-d+1,1) = (u.v(d:end-d+1,2:p+1) * b)/25;

    % Right side
    n_(d:end-d+1,end) = (n(d:end-d+1,end-p:end-1) * b)/25;
    u_.v(d:end-d+1,end) = (u.v(d:end-d+1,end-p:end-1) * b)/25;

    % Bot side
    n_(1,d:end-d+1) = (b' * n(2:p+1,d:end-d+1))/25;
    u_.u(1,d:end-d+1) = (b' *u_.u(2:p+1,d:end-d+1))/25;

    % Top side
    n_(end,d:end-d+1) = (b' * n(end-p:end-1,d:end-d+1))/25;
    u_.u(end,d:end-d+1) = (b' * u.u(end-p:end-1,d:end-d+1))/25;
    

    % Left
    u_.v(d:end-d+1,1) = 0;
    % Right
    u_.v(d:end-d+1,end) = 0;
    % Top
    u_.u(1,d:end-d+1) = 0;
    % Bot
    u_.u(end,d:end-d+1) = 0;
end



















