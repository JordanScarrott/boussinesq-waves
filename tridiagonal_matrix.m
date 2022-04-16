
% Diagonals: a,b,c
% Matrix size: m x m matrix
function tridiag_mat = tridiagonal_matrix(a,b,c)
    dims = size(a);
    tridiag_mat = zeros(dims(2), dims(2));
    
    
    tridiag_mat(1,1) = b(1);
    
    for i=2:dims(2)
        tridiag_mat(i,i-1) = a(i);
        tridiag_mat(i,i) = b(i);
        tridiag_mat(i-1,i) = c(i);
    end
end

