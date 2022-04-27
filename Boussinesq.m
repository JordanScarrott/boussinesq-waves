classdef Boussinesq
    %BOUSSINESQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        animate;
        T;

        iterations;
        tol;
        A0;
        h0;
        h;
        dx;
        dy;
        dt;
        real_x;
        real_y;

        limit;
        filtering;
        filter_period;
        boundary_depth;

        scale;
        xmax;
        ymax;

        x;
        y;
        t;

        xn;
        yn;

        floorProfile;
        initialCondition;

        beta;
        za;
        a1;
        a2;
        b1;
        b2;
        g;

        n;
        nt;
        u;
        ut;
        U;
        V;
        E;
        F;
        F1;
        G;
        G1;

        n_est;
        u_est;
        v_est;
        
        n_error;
        u_error;
        v_error;
        corrector_count;

        u_coeff_mat;
        v_coeff_mat;
    end
    
    methods
        function obj = Boussinesq(setup_args)
            % Number of iterations
            obj.iterations = setup_args(1);
            % Tolerance and max iterations for the corrector step
            obj.tol = setup_args(2);
            obj.limit = 50;
            % Apply filtering
            obj.filtering = 1;
            obj.filter_period = 50;
            % Boundaries
            obj.boundary_depth = 2;


            % Amplitude of the driven waves
            obj.A0 = setup_args(3);
            % Resting height of the water surface
            obj.h0 = setup_args(4);

            % Spacial and Time steps
            obj.dx = setup_args(5);
            obj.dy = setup_args(6);
            obj.dt = setup_args(7);

            % Real length
            obj.real_x = setup_args(8);
            obj.real_y = setup_args(9);


            % Temporary spatial scaling factor
            obj.scale = 1;
            % Max steps in x direction
            obj.xmax = obj.real_x*obj.scale;
            obj.ymax = obj.real_y*obj.scale;
            
            % Discrete Dimensions
            obj.x = 0:obj.dx:obj.xmax;
            obj.y = 0:obj.dy:obj.ymax;
            obj.t = 0:obj.dt:obj.iterations*obj.dt;
            
            % Number of cells in the x and y directions
            obj.xn = numel(obj.x);
            obj.yn = numel(obj.y);

            % Water depth (floor profile)
            obj.floorProfile = FloorProfile(setup_args(10), obj.x, obj.h0);
            obj.h = obj.floorProfile.y_data .* ones(obj.yn, obj.xn);

            % Constants
            obj.beta = -0.531;
            obj.za = obj.beta * obj.h;
            obj.a1 = obj.beta.^2 ./ 2 - 1/6;
            obj.a2 = obj.beta + 1/2;
            obj.b1 = obj.beta.^2 ./ 2;
            obj.b2 = obj.beta;
            obj.g = 9.81;
            
            % Surface Elevation & Horizontal Velocities and Velocity-Potentials
            obj.n = zeros(obj.yn, obj.xn, obj.iterations);
            obj.nt = zeros(obj.yn,obj.xn, obj.iterations);
            obj.u = Vector_Volume([obj.yn, obj.xn], obj.iterations+1);
            obj.ut = Vector_Volume([obj.yn, obj.xn], obj.iterations+1);
            obj.U = zeros(obj.yn, obj.xn, obj.iterations);
            obj.V = zeros(obj.yn, obj.xn, obj.iterations);
            % Potentials
            obj.E = zeros(obj.yn, obj.xn, obj.iterations);
            obj.F = zeros(obj.yn, obj.xn, obj.iterations);
            obj.F1 = zeros(obj.yn, obj.xn, obj.iterations);
            obj.G = zeros(obj.yn, obj.xn, obj.iterations);
            obj.G1 = zeros(obj.yn, obj.xn, obj.iterations);
            
            % Adams-Moulton Estimates: n_est[current, old]
            obj.n_est = zeros(obj.yn, obj.xn, 2);
            obj.u_est = zeros(obj.yn, obj.xn, 2);
            obj.v_est = zeros(obj.yn, obj.xn, 2);
            % Adams-Moulton errors
            obj.n_error = inf * ones(1,obj.iterations);
            obj.u_error = inf * ones(1,obj.iterations);
            obj.v_error = inf * ones(1,obj.iterations);
            obj.corrector_count = zeros(1,obj.iterations);
            
            
            % Precompute coefficient matrices
            obj.u_coeff_mat = u_coeff_matrices(obj.h, obj.b1, obj.b2, obj.dx, obj.xn, obj.yn);
            obj.v_coeff_mat = v_coeff_matrices(obj.h, obj.b1, obj.b2, obj.dy, obj.xn, obj.yn);
            
            % Initial Conditions
            obj.initialCondition = InitialCondition(setup_args(11), obj.n(:,:,1), obj.A0, obj.x, obj.y);
            obj.n(:,:,1) = obj.initialCondition.n;
        end

        function obj = solve(obj)
            for i=1:obj.iterations
                i
                
                if i == 1
                    % Get E, F, G, F1, and G1 for i=1
                    obj.E(:,:,i) = Compute.E(obj.n(:,:,i), obj.u(i).u, obj.u(i).v, obj.h, [obj.a1 obj.a2 obj.b1 obj.b2 obj.g obj.dx obj.dy]);
                    obj.F(:,:,i) = Compute.F(obj.n(:,:,i), obj.u(i).u, obj.u(i).v, obj.g, obj.dx, obj.dy);
                    obj.G(:,:,i) = Compute.G(obj.n(:,:,i), obj.u(i).u, obj.u(i).v, obj.g, obj.dx, obj.dy);
                    obj.F1(:,:,i) = Compute.F1(obj.h, obj.u(i).v, obj.dx, obj.dy, obj.b1, obj.b2);
                    obj.G1(:,:,i) = Compute.G1(obj.h, obj.u(i).u, obj.dx, obj.dy, obj.b1, obj.b2);
            
                    % PREDICTOR - 1st Order AB
                    obj.n(:,:,i+1) = obj.n(:,:,i) + obj.dt * obj.E(:,:,i);
                    obj.U(:,:,i+1) = obj.U(:,:,i) + obj.dt * obj.F(:,:,i);
                    obj.V(:,:,i+1) = obj.V(:,:,i) + obj.dt * obj.G(:,:,i);
                elseif i == 2
                    % PREDICTOR - 2nd Order AB
                    obj.n(:,:,i+1) = obj.n(:,:,i) + obj.dt/2 * (3*obj.E(:,:,i) - obj.E(:,:,i-1));
                    obj.U(:,:,i+1) = obj.U(:,:,i) + obj.dt/2 * (3*obj.F(:,:,i) - obj.F(:,:,i-1)) + obj.F(:,:,i) - obj.F(:,:,i-1);
                    obj.V(:,:,i+1) = obj.V(:,:,i) + obj.dt/2 * (3*obj.G(:,:,i) - obj.G(:,:,i-1)) + obj.G(:,:,i) - obj.G(:,:,i-1);
                else
                    % PREDICTOR SCHEME - 3rd Order AB scheme
                    obj.n(:,:,i+1) = obj.n(:,:,i) + obj.dt/12 * (23*obj.E(:,:,i) - 16*obj.E(:,:,i-1) + 5*obj.E(:,:,i-2));
                    obj.U(:,:,i+1) = obj.U(:,:,i) + obj.dt/12 * (23*obj.F(:,:,i) - 16*obj.F(:,:,i-1) + 5*obj.F(:,:,i-2)) + 2*obj.F1(:,:,i) - 3*obj.F1(:,:,i-1) + obj.F1(:,:,i-2);
                    obj.V(:,:,i+1) = obj.V(:,:,i) + obj.dt/12 * (23*obj.G(:,:,i) - 16*obj.G(:,:,i-1) + 5*obj.G(:,:,i-2)) + 2*obj.G1(:,:,i) - 3*obj.G1(:,:,i-1) + obj.G1(:,:,i-2);
                end
            
                % Compute u and v
                obj.u(i+1).u = solve_for_u(obj.u_coeff_mat, obj.U(:,:,i+1));
                obj.u(i+1).v = solve_for_v(obj.v_coeff_mat, obj.V(:,:,i+1));
            
                % Add boundary conditions
                [obj.n(:,:,i+1), obj.u(i+1)] = boundary_cond(obj.n(:,:,i+1), obj.u(i+1),obj.boundary_depth,obj.dx,obj.dy,obj.A0,obj.t(i),obj.dt);
                
                
                % Check if we need to keep iterating the corrector
                while ((obj.n_error(i) > obj.tol || obj.u_error(i) > obj.tol || obj.v_error(i) > obj.tol))
                    % Track errors for the corrector step
                    obj.n_est(:,:,2) = obj.n(:,:,i+1);
                    obj.u_est(:,:,2) = obj.u(i+1).u;
                    obj.v_est(:,:,2) = obj.u(i+1).v;
                    obj.corrector_count(i) = obj.corrector_count(i) + 1;
            
                    % Compute E, F, G, F1, G1 
                    obj.E(:,:,i+1) = Compute.E(obj.n(:,:,i+1), obj.u(i+1).u, obj.u(i+1).v, obj.h, [obj.a1 obj.a2 obj.b1 obj.b2 obj.g obj.dx obj.dy]);
                    obj.F(:,:,i+1) = Compute.F(obj.n(:,:,i+1), obj.u(i+1).u, obj.u(i+1).v, obj.g, obj.dx, obj.dy);
                    obj.G(:,:,i+1) = Compute.G(obj.n(:,:,i+1), obj.u(i+1).u, obj.u(i+1).v, obj.g, obj.dx, obj.dy);
                    obj.F1(:,:,i+1) = Compute.F1(obj.h, obj.u(i+1).v, obj.dx, obj.dy, obj.b1, obj.b2);
                    obj.G1(:,:,i+1) = Compute.G1(obj.h, obj.u(i+1).u, obj.dx, obj.dy, obj.b1, obj.b2);
            
                    if i == 1
                        % CORRECTOR - 2nd Order AM
                        obj.n(:,:,i+1) = obj.n(:,:,i) + obj.dt/2 * (obj.E(:,:,i+1) + obj.E(:,:,i));
                        obj.U(:,:,i+1) = obj.U(:,:,i) + obj.dt/2 * (obj.F(:,:,i+1) + obj.F(:,:,i)) + obj.F1(:,:,i+1) - obj.F1(:,:,i);
                        obj.V(:,:,i+1) = obj.V(:,:,i) + obj.dt/2 * (obj.G(:,:,i+1) + obj.G(:,:,i)) + obj.G1(:,:,i+1) - obj.G1(:,:,i);
                    elseif i ==2
                        % CORRECTOR - 3rd Order AM
                        obj.n(:,:,i+1) = obj.n(:,:,i) + obj.dt/12 * (5*obj.E(:,:,i+1) + 8*obj.E(:,:,i) - obj.E(:,:,i-1));
                        obj.U(:,:,i+1) = obj.U(:,:,i) + obj.dt/12 * (5*obj.F(:,:,i+1) + 8*obj.F(:,:,i) - obj.F(:,:,i-1)) + obj.F1(:,:,i+1) - obj.F1(:,:,i);
                        obj.V(:,:,i+1) = obj.V(:,:,i) + obj.dt/12 * (5*obj.G(:,:,i+1) + 8*obj.G(:,:,i) - obj.G(:,:,i-1)) + obj.G1(:,:,i+1) - obj.G1(:,:,i);
                    else
                        % CORRECTOR SCHEME - 4th Order AM scheme
                        obj.n(:,:,i+1) = obj.n(:,:,i) + obj.dt/24 * (9*obj.E(:,:,i+1) + 19*obj.E(:,:,i) - 5*obj.E(:,:,i-1) + obj.E(:,:,i-2));
                        obj.U(:,:,i+1) = obj.U(:,:,i) + obj.dt/24 * (9*obj.F(:,:,i+1) + 19*obj.F(:,:,i) - 5*obj.F(:,:,i-1) + obj.F(:,:,i-2)) + obj.F1(:,:,i+1) - obj.F1(:,:,i);
                        obj.V(:,:,i+1) = obj.V(:,:,i) + obj.dt/24 * (9*obj.G(:,:,i+1) + 19*obj.G(:,:,i) - 5*obj.G(:,:,i-1) + obj.G(:,:,i-2)) + obj.G1(:,:,i+1) - obj.G1(:,:,i);
                    end
                    
                    % Compute u and v
                    obj.u(i+1).u = solve_for_u(obj.u_coeff_mat, obj.U(:,:,i+1));
                    obj.u(i+1).v = solve_for_v(obj.v_coeff_mat, obj.V(:,:,i+1));
                    
                    % Add boundary conditions
                    [obj.n(:,:,i+1), obj.u(i+1)] = boundary_cond(obj.n(:,:,i+1), obj.u(i+1),obj.boundary_depth,obj.dx,obj.dy,obj.A0,obj.t(i),obj.dt);

                    % Store estimates for this iteration so we can compute error
                    obj.n_est(:,:,1) = obj.n(:,:,i+1);
                    obj.u_est(:,:,1) = obj.u(i+1).u;
                    obj.v_est(:,:,1) = obj.u(i+1).v;
            
                    % Compute error for n, u, and v
                    obj.n_error(i) = sum(abs(obj.n_est(:,:,1) - obj.n_est(:,:,2)), [1 2]) / sum(abs(obj.n_est(:,:,1)), [1 2]);
                    obj.u_error(i) = sum(abs(obj.u_est(:,:,1) - obj.u_est(:,:,2)), [1 2]) / sum(abs(obj.u_est(:,:,1)), [1 2]);
                    obj.v_error(i) = sum(abs(obj.v_est(:,:,1) - obj.v_est(:,:,2)), [1 2]) / sum(abs(obj.v_est(:,:,1)), [1 2]);
                    
                    clc
                    fprintf('corrector_count = %.d\n', obj.corrector_count(i))
                    fprintf('n_error = %.10f\n', obj.n_error(i))
                    fprintf('u_error = %.10f\n', obj.u_error(i))
                    fprintf('v_error = %.10f\n', obj.v_error(i))
                    
                    if (obj.corrector_count(i) == obj.limit)
                        fprintf('\n\nProgram terminated at:\nIteration: %d \nCorrector step: %d\n', i, obj.corrector_count(i))
                        fprintf('Errors this iteration:\nn_error: %.10f \nu_error: %.10f \nv_error: %.10f\n', obj.n_error(i), obj.u_error(i), obj.v_error(i))
                        error('Non-convergeance')
                    end
                end
                
                % Update E, F, G, F1, and G1 for the finalized i+1
                obj.E(:,:,i+1) = Compute.E(obj.n(:,:,i+1), obj.u(i+1).u, obj.u(i+1).v, obj.h, [obj.a1 obj.a2 obj.b1 obj.b2 obj.g obj.dx obj.dy]);
                obj.F(:,:,i+1) = Compute.F(obj.n(:,:,i+1), obj.u(i+1).u, obj.u(i+1).v, obj.g, obj.dx, obj.dy);
                obj.G(:,:,i+1) = Compute.G(obj.n(:,:,i+1), obj.u(i+1).u, obj.u(i+1).v, obj.g, obj.dx, obj.dy);
                obj.F1(:,:,i+1) = Compute.F1(obj.h, obj.u(i+1).v, obj.dx, obj.dy, obj.b1, obj.b2);
                obj.G1(:,:,i+1) = Compute.G1(obj.h, obj.u(i+1).u, obj.dx, obj.dy, obj.b1, obj.b2);
            
                % Move the current estimates to the old estimate slots
                obj.n_est(:,:,2) = obj.n_est(:,:,1);
                obj.u_est(:,:,2) = obj.u_est(:,:,1);
                obj.v_est(:,:,2) = obj.v_est(:,:,1);  

                if (obj.filtering == 1)
                    if rem(i,obj.filter_period) == 0
                        obj.n(:,:,i+1) = filter2d(obj.n(:,:,i+1));
                        obj.u(i+1).u = filter2d(obj.u(i+1).u);
                        obj.u(i+1).v = filter2d(obj.u(i+1).v);
                    end
                end
            end
        end

        % Overload function for displaySpecificMeshes
        function obj = displayMeshes(obj)
            obj = obj.displaySpecificMeshes(1);
        end

        function obj = displaySpecificMeshes(obj, percentageOfTotalIterations)
            iterationsToDisplay = round(percentageOfTotalIterations * obj.iterations);

            obj.animate = 1;
            obj.T = 5;

            % Plotting Final Meshes
            % chart_titles = ["n(:,:,i)", "U(:,:,i)", "V(:,:,i)", "E(:,:,i)", "F(:,:,i)", "G(:,:,i)"];
%             chart_titles = ["n(:,:,i)", "U(:,:,i)", "u(:).u", "E(:,:,i)", "F(:,:,i)", "u(:).v"];
            chart_titles = ["n(:,:,i)"];
            % display_meshes(obj.animate, obj.T, cat(3,obj.n,obj.U,obj.V,obj.E,obj.F,obj.G), obj.iterations+1, chart_titles)
            % display_meshes(obj.animate, obj.T, cat(3,obj.n,obj.U,obj.V,obj.E), obj.iterations+1, chart_titles)
%             display_meshes(obj.animate, obj.T, cat(3, obj.n, obj.U, obj.u(:).u, obj.E, obj.F, obj.u(:).v), obj.iterations+1, chart_titles)
            display_meshes(obj.animate, obj.T, obj.n, iterationsToDisplay+1, chart_titles);


            % Plotting in one dimension
            for i=1:round(iterationsToDisplay/5):iterationsToDisplay
                figure(10)
                plot(obj.n(round(obj.yn/2),:,i))
                ylim([-0.05 0.06])
                pause(0.1)
            end
            plot(obj.n(round(obj.yn/2),:,end))
            ylim([-0.05 0.06])

            figure(3)
            contour(obj.n(:,:,iterationsToDisplay+1))
        end

        function obj = saveParamData(obj)
            filter = {'.xlsx', 'Excel Sheet (.xlsx)';'*.*', 'All Files (*.*)'};
            [file,path,index] = uiputfile(filter, "Save File Name");
            
            dataNames = ["iterations", "tol", "A0", "h0", "dx", "dy", "dt", "real_x", "real_y", "limit", "filtering", "filter_period", "boundary_depth", "scale", "xmax", "ymax", "xn", "yn", "beta", "a1", "a2", "b1", "b2", "g", "obj.floorProfile.SELECTION", "obj.initialCondition.SELECTION"];
            dataToSave = [obj.iterations, obj.tol, obj.A0, obj.h0, obj.dx, obj.dy, obj.dt, obj.real_x, obj.real_y, obj.limit, obj.filtering, obj.filter_period, obj.boundary_depth, obj.scale, obj.xmax, obj.ymax, obj.xn, obj.yn, obj.beta, obj.a1, obj.a2, obj.b1, obj.b2, obj.g];
            textData = [obj.floorProfile.SELECTION, obj.initialCondition.SELECTION];

            dataToSave = [dataNames; dataToSave textData];
            writematrix(dataToSave, strcat(path, file))
        end

    end

    % Here are the static methods for the Boussinesq class
    methods(Static)
        function A = loadPresetFromFile(obj)
            filter = {'.xlsx', 'Excel Sheet (.xlsx)';'*.*', 'All Files (*.*)'};
            [file,path,index] = uigetfile(filter);

            tableData = table2cell(readtable(strcat(path, file)));

            numberData = str2double(tableData(1,1:9));
            textData = string(tableData(1,end-1:end));

            functionParams = [numberData, FloorProfile.getProfileEnum(textData(1,1)), InitialCondition.getProfileEnum(textData(1,2))];

            A = Boussinesq(functionParams);
        end

    end
end

