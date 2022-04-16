function ui = solve_for_u(u_coeff_mat,Ui)
    dims = size(Ui);
    ui = zeros(dims(1),dims(2));
        
%     if rank(Ui) > 0 
        % For each row, k
        for k=1:dims(1)
            ui(k,:) = tridiagonal(u_coeff_mat(:,:,k), Ui(k,:));
        end
%     end
    

%     % Test function
%     soln  = zeros(dims(1),dims(2));
%     for k=1:dims(1)
%         soln(k,:) = u_coeff_mat(:,:,k) * ui(k,:)';
%     end

    
end

