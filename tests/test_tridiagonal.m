% Test for tridiagonal.m

% Define a tridiagonal matrix
A = [2, -1, 0;
     -1, 2, -1;
     0, -1, 2];

% Define the vector d
d = [1; 2; 3];

% Calculate the expected solution
expected_x = A\d;

% Call the function
actual_x = tridiagonal(A, d);

% Compare the actual and expected solutions
assert(isequal(actual_x, expected_x), 'Test for tridiagonal.m failed');

fprintf('Test for tridiagonal.m passed.\n');
