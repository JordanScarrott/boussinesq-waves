classdef Profile_description
    %PROFILE_DESCRIPTION
    %   Generates a number of points on a line between 0 and 1
    
    properties
        % The 2D point array and x and y components
        points;
        x_points;
        y_points;

        % The new higher resolution x and y arrays
        x;
        y;

        % All x related quantities
        dx;
        samples;
        segments;
        num_of_segments;
    end
    
    methods
        function obj = Profile_description(points_, x_)
            obj.points = points_;
            obj.x_points = obj.points(:,1);
            obj.y_points = obj.points(:,2);

            obj.x = x_;
            obj.y = zeros(1, length(obj.x));

            obj.dx = obj.x(2) - obj.x(1);
            obj.samples = length(obj.x)-1;
            obj.segments = floor(obj.x_points .* obj.samples) + 1;
            obj.num_of_segments = length(obj.x_points)-1;

            % For each segment fill in the missing data
            for i=1:obj.num_of_segments
                % Delta y for current segment
                delta_y = obj.y_points(i+1) - obj.y_points(i);
                num_of_x_intervals = obj.segments(i+1) - obj.segments(i);
                
                grad_per_dx_for_segment = delta_y / (num_of_x_intervals * obj.dx);

                intervals = (1:num_of_x_intervals);
                    
                obj.y(obj.segments(i)+1:obj.segments(i+1)) = obj.y(obj.segments(i)) + grad_per_dx_for_segment .* intervals * obj.dx;
            end

            % This converted the floor height above zero into water depth
            obj.y = (obj.y - 1) * -1;
        end
    end
end

