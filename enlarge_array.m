
% Simply adds white space to the edges of an array to make it a certain
% size. ONLY WORKS FOR ENLARGING ARRAYS NOT SHRINKING (rows, cols) = (m, n)
function enlarged_array = enlarge_array(A, m, n)
    dims = size(A);
    enlarged_array = A;
    
    if numel(dims) == 3
        % If it is a 3D array
        enlarged_array(dims(1)+1:m,dims(2)+1:n,:) = zeros(m-dims(1), n-dims(2), dims(3));
    elseif numel(dims) == 2
        % Else if its a 2D array
        enlarged_array(dims(1)+1:m,dims(2)+1:n) = zeros(m-dims(1), n-dims(2));
    elseif numel(dims) == 1
    % Else if 1D array
        enlarged_array(dims(2)+1:n) = zeros(1, n-dims(2));
    end
end

