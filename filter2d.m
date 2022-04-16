function filtered_data = filter2d(n)
    
    dims = size(n);
    filtered_data = n;
    temp = filtered_data;
    

    % First in the x direction for all y values
    for k=5:dims(2)-4
        temp(:,k) = 1/256 * (186 * n(:,k) + 56 * (n(:,k+1) + n(:,k-1)) - 28 * (n(:,k+2)+ n(:,k-2)) + 8 * (n(:,k+3) + n(:,k-3)) - (n(:,k+4) + n(:,k-4)));
    end
    
    % Then in the y direction for all x values
    for j=5:dims(1)-4
        filtered_data(j,:) = 1/256 * (186 * temp(j,:) + 56 * (temp(j+1,:) + temp(j-1,:)) - 28 * (temp(j+2,:)+ temp(j-2,:)) + 8 * (temp(j+3,:) + temp(j-3,:)) - (temp(j+4,:) + temp(j-4,:)));
    end
    
    
    % And now a 3rd Order Filter for the Boundaries
    % First in the x direction for all y values
    k=4;
    temp(:,k) = 1/64 * (44 * n(:,k) + 15 * (n(:,k+1) + n(:,k-1)) - 6 * (n(:,k+2)+ n(:,k-2)) + 1 * (n(:,k+3) + n(:,k-3)));
    k=dims(2)-3;
    temp(:,k) = 1/64 * (44 * n(:,k) + 15 * (n(:,k+1) + n(:,k-1)) - 6 * (n(:,k+2)+ n(:,k-2)) + 1 * (n(:,k+3) + n(:,k-3)));
    
    % Then in the y direction for all x values
    j=4;
    filtered_data(j,:) = 1/64 * (44 * temp(j,:) + 15 * (temp(j+1,:) + temp(j-1,:)) - 6 * (temp(j+2,:)+ temp(j-2,:)) + 1 * (temp(j+3,:) + temp(j-3,:)));
    j=dims(1)-3;
    filtered_data(j,:) = 1/64 * (44 * temp(j,:) + 15 * (temp(j+1,:) + temp(j-1,:)) - 6 * (temp(j+2,:)+ temp(j-2,:)) + 1 * (temp(j+3,:) + temp(j-3,:)));
    
    
    % And now a 2nd Order Filter for the Boundaries
    % First in the x direction for all y values
    k=3;
    temp(:,k) = 1/16 * (10 * n(:,k) + 4 * (n(:,k+1) + n(:,k-1)) - 1 * (n(:,k+2)+ n(:,k-2)));
    k=dims(2)-2;
    temp(:,k) = 1/16 * (10 * n(:,k) + 4 * (n(:,k+1) + n(:,k-1)) - 1 * (n(:,k+2)+ n(:,k-2)));
    
    % Then in the y direction for all x values
    j=3;
    filtered_data(j,:) = 1/16 * (10 * temp(j,:) + 4 * (temp(j+1,:) + temp(j-1,:)) - 1 * (temp(j+2,:)+ temp(j-2,:)));
    j=dims(1)-2;
    filtered_data(j,:) = 1/16 * (10 * temp(j,:) + 4 * (temp(j+1,:) + temp(j-1,:)) - 1 * (temp(j+2,:)+ temp(j-2,:)));
    
    
    % And now a 1st Order Filter for the Boundaries
    % First in the x direction for all y values
    k=2;
    temp(:,k) = 1/4 * (2 * n(:,k) + 1 * (n(:,k+1) + n(:,k-1)));
    k=dims(2)-1;
    temp(:,k) = 1/4 * (2 * n(:,k) + 1 * (n(:,k+1) + n(:,k-1)));
    
    % Then in the y direction for all x values
    j=2;
    filtered_data(j,:) = 1/4 * (2 * temp(j,:) + 1 * (temp(j+1,:) + temp(j-1,:)));
    j=dims(1)-1;
    filtered_data(j,:) = 1/4 * (2 * temp(j,:) + 1 * (temp(j+1,:) + temp(j-1,:)));
    
    
    
end

