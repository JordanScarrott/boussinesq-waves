
% Jordan Scarrott
% 215179641
% 15 April 2022

% Adams Bashfourth Predictor Corrector method for the Boussinesq Equation
% Wei 1995

% Changelog:
%   - This version focuses on abstracting this file into a function that
%   can be run to solve the Boussinesq equation more easily

clear
clc
animate = 1;
T = 5;

% Number of iterations
iterations = 100;
% Tolerance and max iterations for the corrector step
tol = 0.001;
limit = 50;
% Apply filtering
filtering = 1;
filter_period = 50;
% Boundaries
boundary_depth = 2;

% Amplitude of the driven waves
A0 = 0.045;
% Resting height of the water surface
h0 = 0.45;

% Time step
dx = 0.075;
dy = 0.075;
dt = 0.05;

% Real length
real_x = 7.5;
real_y = 7.5;
% Temporary spatial scaling factor
scale = 1;
% Max steps in x direction
xmax = real_x*scale;
ymax = real_y*scale;

% Discrete Dimensions
x = 0:dx:xmax;
y = 0:dy:ymax;
t = 0:dt:iterations*dt;

% Number of cells in the x and y directions
xn = numel(x);
yn = numel(y);

% The water depth, h
% h = h0-h0/(real_x*12.7/13) * ((ones(yn, xn) .* x ./ scale) - (real_x*0.3/13));
h = h0 * ones(yn, xn);
% floor_profile = [[0, 0]; [0.25, 0]; [0.35, h0/3]; [0.4, h0/3]; [0.45, h0/5]; [0.5, h0/5]; [0.75, h0]; [1, h0]];
% h = Profile_description(floor_profile, x).y .* ones(yn,xn);

% Constants
beta = -0.531;
za = beta * h;
a1 = beta.^2 ./ 2 - 1/6;
a2 = beta + 1/2;
b1 = beta.^2 ./ 2;
b2 = beta;
g = 9.81;

% Surface Elevation & Horizontal Velocities and Velocity-Potentials
n = zeros(yn, xn, iterations);
nt = zeros(yn,xn, iterations);
u = Vector_Volume([yn, xn], iterations+1);
ut = Vector_Volume([yn, xn], iterations+1);
U = zeros(yn, xn, iterations);
V = zeros(yn, xn, iterations);
% Potentials
E = zeros(yn, xn, iterations);
F = zeros(yn, xn, iterations);
F1 = zeros(yn, xn, iterations);
G = zeros(yn, xn, iterations);
G1 = zeros(yn, xn, iterations);

% Adams-Moulton Estimates: n_est[current, old]
n_est = zeros(yn, xn, 2);
u_est = zeros(yn, xn, 2);
v_est = zeros(yn, xn, 2);
% Adams-Moulton errors
n_error = inf * ones(1,iterations);
u_error = inf * ones(1,iterations);
v_error = inf * ones(1,iterations);
corrector_count = zeros(1,iterations);


% Precompute coefficient matrices
u_coeff_mat = u_coeff_matrices(h,b1,b2,dx,xn,yn);
v_coeff_mat = v_coeff_matrices(h,b1,b2,dy,xn,yn);

% Initial Conditions
n(:,:,1) = initial_condition(n(:,:,1),A0,x,y);
% u(1).u = initial_sech(u(1).u,A0,1);

for i=1:iterations
    i
    
    if i == 1
        % Get E, F, G, F1, and G1 for i=1
        E(:,:,i) = Compute.E(n(:,:,i), u(i).u, u(i).v, h, [a1 a2 b1 b2 g dx dy]);
        F(:,:,i) = Compute.F(n(:,:,i), u(i).u, u(i).v, g, dx, dy);
        G(:,:,i) = Compute.G(n(:,:,i), u(i).u, u(i).v, g, dx, dy);
        F1(:,:,i) = Compute.F1(h, u(i).v, dx, dy, b1, b2);
        G1(:,:,i) = Compute.G1(h, u(i).u, dx, dy, b1, b2);

        % PREDICTOR - 1st Order AB
        n(:,:,i+1) = n(:,:,i) + dt * E(:,:,i);
        U(:,:,i+1) = U(:,:,i) + dt * F(:,:,i);
        V(:,:,i+1) = V(:,:,i) + dt * G(:,:,i);
    elseif i == 2
        % PREDICTOR - 2nd Order AB
        n(:,:,i+1) = n(:,:,i) + dt/2 * (3*E(:,:,i) - E(:,:,i-1));
        U(:,:,i+1) = U(:,:,i) + dt/2 * (3*F(:,:,i) - F(:,:,i-1)) + F(:,:,i) - F(:,:,i-1);
        V(:,:,i+1) = V(:,:,i) + dt/2 * (3*G(:,:,i) - G(:,:,i-1)) + G(:,:,i) - G(:,:,i-1);
    else
        % PREDICTOR SCHEME - 3rd Order AB scheme
        n(:,:,i+1) = n(:,:,i) + dt/12 * (23*E(:,:,i) - 16*E(:,:,i-1) + 5*E(:,:,i-2));
        U(:,:,i+1) = U(:,:,i) + dt/12 * (23*F(:,:,i) - 16*F(:,:,i-1) + 5*F(:,:,i-2)) + 2*F1(:,:,i) - 3*F1(:,:,i-1) + F1(:,:,i-2);
        V(:,:,i+1) = V(:,:,i) + dt/12 * (23*G(:,:,i) - 16*G(:,:,i-1) + 5*G(:,:,i-2)) + 2*G1(:,:,i) - 3*G1(:,:,i-1) + G1(:,:,i-2);
    end

    % Compute u and v
    u(i+1).u = solve_for_u(u_coeff_mat, U(:,:,i+1));
    u(i+1).v = solve_for_v(v_coeff_mat, V(:,:,i+1));

    % Add boundary conditions
    [n(:,:,i+1), u(i+1)] = boundary_cond(n(:,:,i+1), u(i+1),boundary_depth,dx,dy,A0,t(i),dt);
    
    
    % Check if we need to keep iterating the corrector
    while ((n_error(i) > tol || u_error(i) > tol || v_error(i) > tol))
        % Track errors for the corrector step
        n_est(:,:,2) = n(:,:,i+1);
        u_est(:,:,2) = u(i+1).u;
        v_est(:,:,2) = u(i+1).v;
        corrector_count(i) = corrector_count(i) + 1;

        % Compute E, F, G, F1, G1 
        E(:,:,i+1) = Compute.E(n(:,:,i+1), u(i+1).u, u(i+1).v, h, [a1 a2 b1 b2 g dx dy]);
        F(:,:,i+1) = Compute.F(n(:,:,i+1), u(i+1).u, u(i+1).v, g, dx, dy);
        G(:,:,i+1) = Compute.G(n(:,:,i+1), u(i+1).u, u(i+1).v, g, dx, dy);
        F1(:,:,i+1) = Compute.F1(h, u(i+1).v, dx, dy, b1, b2);
        G1(:,:,i+1) = Compute.G1(h, u(i+1).u, dx, dy, b1, b2);

        if i == 1
            % CORRECTOR - 2nd Order AM
            n(:,:,i+1) = n(:,:,i) + dt/2 * (E(:,:,i+1) + E(:,:,i));
            U(:,:,i+1) = U(:,:,i) + dt/2 * (F(:,:,i+1) + F(:,:,i)) + F1(:,:,i+1) - F1(:,:,i);
            V(:,:,i+1) = V(:,:,i) + dt/2 * (G(:,:,i+1) + G(:,:,i)) + G1(:,:,i+1) - G1(:,:,i);
        elseif i ==2
            % CORRECTOR - 3rd Order AM
            n(:,:,i+1) = n(:,:,i) + dt/12 * (5*E(:,:,i+1) + 8*E(:,:,i) - E(:,:,i-1));
            U(:,:,i+1) = U(:,:,i) + dt/12 * (5*F(:,:,i+1) + 8*F(:,:,i) - F(:,:,i-1)) + F1(:,:,i+1) - F1(:,:,i);
            V(:,:,i+1) = V(:,:,i) + dt/12 * (5*G(:,:,i+1) + 8*G(:,:,i) - G(:,:,i-1)) + G1(:,:,i+1) - G1(:,:,i);
        else
            % CORRECTOR SCHEME - 4th Order AM scheme
            n(:,:,i+1) = n(:,:,i) + dt/24 * (9*E(:,:,i+1) + 19*E(:,:,i) - 5*E(:,:,i-1) + E(:,:,i-2));
            U(:,:,i+1) = U(:,:,i) + dt/24 * (9*F(:,:,i+1) + 19*F(:,:,i) - 5*F(:,:,i-1) + F(:,:,i-2)) + F1(:,:,i+1) - F1(:,:,i);
            V(:,:,i+1) = V(:,:,i) + dt/24 * (9*G(:,:,i+1) + 19*G(:,:,i) - 5*G(:,:,i-1) + G(:,:,i-2)) + G1(:,:,i+1) - G1(:,:,i);
        end
        
        % Compute u and v
        u(i+1).u = solve_for_u(u_coeff_mat, U(:,:,i+1));
        u(i+1).v = solve_for_v(v_coeff_mat, V(:,:,i+1));
        
        % Add boundary conditions
        [n(:,:,i+1), u(i+1)] = boundary_cond(n(:,:,i+1), u(i+1),boundary_depth,dx,dy,A0,t(i),dt);
    
        % Store estimates for this iteration so we can compute error
        n_est(:,:,1) = n(:,:,i+1);
        u_est(:,:,1) = u(i+1).u;
        v_est(:,:,1) = u(i+1).v;

        % Compute error for n, u, and v
        n_error(i) = sum(abs(n_est(:,:,1) - n_est(:,:,2)), [1 2]) / sum(abs(n_est(:,:,1)), [1 2]);
        u_error(i) = sum(abs(u_est(:,:,1) - u_est(:,:,2)), [1 2]) / sum(abs(u_est(:,:,1)), [1 2]);
        v_error(i) = sum(abs(v_est(:,:,1) - v_est(:,:,2)), [1 2]) / sum(abs(v_est(:,:,1)), [1 2]);
        
        clc
        fprintf('corrector_count = %.d\n', corrector_count(i))
        fprintf('n_error = %.10f\n', n_error(i))
        fprintf('u_error = %.10f\n', u_error(i))
        fprintf('v_error = %.10f\n', v_error(i))
        
        if (corrector_count(i) == limit)
            fprintf('\n\nProgram terminated at:\nIteration: %d \nCorrector step: %d\n', i, corrector_count(i))
            fprintf('Errors this iteration:\nn_error: %.10f \nu_error: %.10f \nv_error: %.10f\n', n_error(i), u_error(i), v_error(i))
            error('Non-convergeance')
        end
    end
    
    % Update E, F, G, F1, and G1 for the finalized i+1
    E(:,:,i+1) = Compute.E(n(:,:,i+1), u(i+1).u, u(i+1).v, h, [a1 a2 b1 b2 g dx dy]);
    F(:,:,i+1) = Compute.F(n(:,:,i+1), u(i+1).u, u(i+1).v, g, dx, dy);
    G(:,:,i+1) = Compute.G(n(:,:,i+1), u(i+1).u, u(i+1).v, g, dx, dy);
    F1(:,:,i+1) = Compute.F1(h, u(i+1).v, dx, dy, b1, b2);
    G1(:,:,i+1) = Compute.G1(h, u(i+1).u, dx, dy, b1, b2);

    % Move the current estimates to the old estimate slots
    n_est(:,:,2) = n_est(:,:,1);
    u_est(:,:,2) = u_est(:,:,1);
    v_est(:,:,2) = v_est(:,:,1);  
        
    if (filtering == 1)
        if rem(i,filter_period) == 0
            n(:,:,i+1) = filter2d(n(:,:,i+1));
            u(i+1).u = filter2d(u(i+1).u);
            u(i+1).v = filter2d(u(i+1).v);
        end
    end
    
%     [n(:,:,i+1), u(i+1)] = wavemaker_boundary(n(:,:,i+1), u(i+1),boundary_depth,dx,dy,A0,t(i),dt);
    
end



% Plotting Final Meshes
% chart_titles = ["n(:,:,i)", "U(:,:,i)", "V(:,:,i)", "E(:,:,i)", "F(:,:,i)", "G(:,:,i)"];
% chart_titles = ["n(:,:,i)", "U(:,:,i)", "u(:).u", "E(:,:,i)", "F(:,:,i)", "u(:).v"];
chart_titles = ["n(:,:,i)"];
% display_meshes(animate, T, cat(3,n,U,V,E,F,G), iterations+1, chart_titles)
% display_meshes(animate, T, cat(3,n,U,V,E), iterations+1, chart_titles)
% display_meshes(animate, T, cat(3,n,U,u(:).u,E,F,u(:).v), iterations+1, chart_titles)
display_meshes(animate, T, n, iterations+1, chart_titles)


% Plotting in one dimension
for i=1:round(iterations/5):iterations
    figure(10)
    plot(n(round(yn/2),:,i))
    ylim([-0.05 0.06])
    pause(0.1)
end
plot(n(round(yn/2),:,end))
ylim([-0.05 0.06])

figure(3)
contour(n(:,:,end))





























% comment...






