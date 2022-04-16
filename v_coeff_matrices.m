function v_coeff_mat = v_coeff_matrices(h,b1,b2,dy,xn,yn)

    % Creating the coefficient matrix for solving for v from V
    coeff_vA = 1/dy^2 * h .* (b1 * h + b2 * [h(2,:); h(1:yn-1,:)]);
    coeff_vB = 1 - 2*h.^2*(b1 + b2)/dy^2;
    coeff_vC = 1/dy^2 * h .* (b1 * h + b2 * [h(2:yn,:); h(yn-1,:)]);

    % Create a triadiagonal coefficient matrix for each column, k
    v_coeff_mat = zeros(yn,yn,xn);
    % For col k
    for k=1:xn
        % Setup the tridiagonal matrix of coefficients
        %   B C 0 0 0 0
        %   A B C 0 0 0 
        %   0 A B C 0 0
        %   0 0 A B C 0
        %   0 0 0 A B C
        %   0 0 0 0 A B
        v_coeff_mat(1,1,k) = coeff_vB(1,k);
        for j=2:yn
            i=j;
            v_coeff_mat(j,i-1,k) = coeff_vA(i,k);
            v_coeff_mat(j,i,k) = coeff_vB(i,k);
            v_coeff_mat(j-1,i,k) = coeff_vC(i,k);
        end
    end
end

