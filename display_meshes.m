
% Animates all given meshes over a given time interval with a update 
% frequency thats determined automatically
function display_meshes(animate, T, meshes, iterations, chart_titles)
    % Displays a full square grid when drawing rectangular datasets
    full_grid = 0;
    % The size of each mesh and how many iterations through time they has
    mesh_dims = size(meshes);

    % Determine number of meshes to be drawn
    num_meshes = mesh_dims(3) / iterations;
    
    if (full_grid == 1)
        % Adding in blank space for rectangular grids
        if mesh_dims(2) ~= mesh_dims(1)

            diff = abs(mesh_dims(2) - mesh_dims(1));

            % If the mesh is deeper than it is wide
            if mesh_dims(2) > mesh_dims(1)
                temp_meshes = zeros(mesh_dims(2),mesh_dims(2),mesh_dims(3));

                for i=1:mesh_dims(3)
                    temp_meshes(:,:,i) = [zeros(floor(diff/2),mesh_dims(2)); meshes(:,:,i); zeros(ceil(diff/2),mesh_dims(2))];
                end
            else
                temp_meshes = zeros(mesh_dims(1),mesh_dims(1),mesh_dims(3));

                for i=1:mesh_dims(3)
                    temp_meshes(:,:,i) = [zeros(mesh_dims(1),floor(diff/2)), meshes(:,:,i), zeros(mesh_dims(1),ceil(diff/2))];
                end
            end

            meshes = temp_meshes;
            clear temp_meshes
        end
    end
    
    % Plotting
    if animate==1    
        % Only draw a few of the iterations
        steps = T*5;
        pause_time = T/steps;

        draw_every = iterations/steps;
                
        % Plot each mesh every iteration
        for i=0:steps-1
            % Time the graphing events so we can adjust the pause time
            start = tic;
            
            for j=1:num_meshes
                figure(j)
%                 mesh(meshes(:,:,round(iterations*(j-1) + draw_every*i)+1));
                meshc(meshes(:,:,round(iterations*(j-1) + draw_every*i)+1));
                % set(gca,'View',[-10 10])
%                 zlim([-0.05 0.06])
                title(chart_titles(j))
                ylabel('y')
                xlabel('x')
                zlabel('Amplitude')
            end
            
            fin = toc(start);
            % Display over a constant time interval of T seconds
            if fin < pause_time
                pause(pause_time - fin)
            end
        end
        
        % Finally, draw the last iteration incase we stopped ealy due to rounding
        % Do this by setting animate to 0 and letting the algorithm plot
        % the final iteration there
        animate = 0;
    end
    
    
    if animate==0
        for i=1:num_meshes
            figure(i)
%             mesh(meshes(:,:,iterations*(i)))
            meshc(meshes(:,:,iterations*(i)));
%             zlim([-0.05 0.06])
            title(chart_titles(i))
            ylabel('y')
            xlabel('x')
            zlabel('Amplitude')
        end
    end
end


