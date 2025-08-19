% Test for tridiagonal_matrix.m

% Define the diagonals
a = [2, 3, 4];
b = [1, 5, 9];
c = [6, 7, 8];

% Expected output
expected_matrix = [1, 6, 0;
                   2, 5, 7;
                   0, 3, 9];

% Call the function
actual_matrix = tridiagonal_matrix(a, b, c);

% Compare the actual and expected matrices
assert(isequal(actual_matrix, expected_matrix), 'Test for tridiagonal_matrix.m failed');

fprintf('Test for tridiagonal_matrix.m passed.\n');
