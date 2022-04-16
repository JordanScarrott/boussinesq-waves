
% Takes a 2d array and returns the double forward derivative in either the x
% or y direction. 1 is x direction, 2 is y direction
function derivative = doublederiv2d(u_data, direction, delta)
    
    if direction == 1
        derivative = (u_data(:, 3:end) - 2*u_data(:, 2:end-1) + u_data(:, 1:end-2)) / delta^2;
        % Just set the border derivatives to the same as their neighbour
%         derivative = [derivative(:,1) derivative derivative(:,end)];
        derivative = [derivative(:,1)-(derivative(:,2)-derivative(:,1)) derivative derivative(:,end)+(derivative(:,end)-derivative(:,end-1))];
    elseif direction == 2
        derivative = (u_data(3:end, :) - 2*u_data(2:end-1, :) + u_data(1:end-2, :)) / delta^2;
        % Reclaim correct dimensions
%         derivative = [derivative(1,:); derivative; derivative(end,:)];
        derivative = [derivative(1,:)-(derivative(2,:)-derivative(1,:)); derivative; derivative(end,:)+(derivative(end,:)-derivative(end-1,:))];
    end
end

