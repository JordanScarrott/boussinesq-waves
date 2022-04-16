function u_coeff_mat = u_coeff_matrices(h,b1,b2,dx,xn,yn)
    
    % Creating the coefficient matrix for solving for u from U
    coeff_uA = 1/dx^2 * h .* (b1 * h + b2 * [h(:,2) h(:,1:xn-1)]);
    coeff_uB = 1 - 2*h.^2*(b1 + b2)/dx^2;
    coeff_uC = 1/dx^2 * h .* (b1 * h + b2 * [h(:,2:xn) h(:,xn-1)]);

    % Create a triadiagonal coefficient matrix for each row, k
    u_coeff_mat = zeros(xn,xn,yn);
    % For row k
    for k=1:yn
        % Setup the tridiagonal matrix of coefficients
        %   B C 0 0 0 0
        %   A B C 0 0 0 
        %   0 A B C 0 0
        %   0 0 A B C 0
        %   0 0 0 A B C
        %   0 0 0 0 A B
        u_coeff_mat(1,1,k) = coeff_uB(k,1);
        for j=2:xn
            i=j;
            u_coeff_mat(j,i-1,k) = coeff_uA(k,i);
            u_coeff_mat(j,i,k) = coeff_uB(k,i);
            u_coeff_mat(j-1,i,k) = coeff_uC(k,i);
        end
    end

end

