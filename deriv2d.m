
% Takes a 2d array and returns the forward derivative in either the x
% or y direction. 1 is x direction, 2 is y direction
function derivative = deriv2d(u_data, direction, delta)
    
    if direction == 1
        derivative = (u_data(:, 3:end) - u_data(:, 1:end-2)) / (2*delta);
        % Adding left and right boundary 1st order deriv
        derivative = [(u_data(:,2)-u_data(:,1))/delta, derivative, (u_data(:,end)-u_data(:,end-1))/delta];
    elseif direction == 2
        derivative = (u_data(3:end, :) - u_data(1:end-2, :)) / (2*delta);
        % Adding top and bot boundary 1st order deriv
        derivative = [(u_data(2,:)-u_data(1,:))/delta; derivative; (u_data(end,:)-u_data(end-1,:))/delta];
    end
    
end

