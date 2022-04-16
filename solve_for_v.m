function vi = solve_for_v(v_coeff_mat,Vi)
    dims = size(Vi);
    vi = zeros(dims(1),dims(2));
        
%     if rank(Vi) > 0 
        % For each column, k
        for k=1:dims(2)
%             vi(:,k) = Vi(:,k) \ v_coeff_mat(:,:,k);
            vi(:,k) = tridiagonal(v_coeff_mat(:,:,k), Vi(:,k));
        end
%     end
    
end

