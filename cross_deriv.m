



function cross_derivative = cross_deriv(u_data, dx, dy)
    
    dims = size(u_data);
    cross_derivative = zeros(dims(1),dims(2));

    for j=3:dims(1)-2
        for k=3:dims(2)-2
            cross_derivative(j,k) = (u_data(j+1,k+1) + u_data(j-1,k-1) - u_data(j+1,k-1) - u_data(j-1,k+1)) / (4*dx*dy);
        end
    end
    
end


























% Comment
