# Testing Guide

This document provides a guide for testing the Boussinesq Wave Simulation.

## Running Tests

The tests are located in the `tests` directory. Each test script is a MATLAB/Octave script that can be run independently. To run a test, you need to have MATLAB or Octave installed.

To run a specific test, navigate to the `tests` directory and run the script from the MATLAB or Octave command line. For example, to run the test for `tridiagonal_matrix.m`, you would run the following command:

```
run('test_tridiagonal_matrix.m')
```

Alternatively, you can run all the tests by creating a master script that calls all the individual test scripts.

**Note:** At the time of writing, the testing environment did not have MATLAB or Octave installed, so the tests could not be run to verify their correctness.

## Adding New Tests

To add a new test, create a new file in the `tests` directory. The file should be named `test_<function_name>.m`, where `<function_name>` is the name of the function you are testing.

The test script should:

1.  Set up the necessary input data for the function being tested.
2.  Call the function with the input data.
3.  Compare the actual output with the expected output.
4.  Use the `assert` function to check if the actual output matches the expected output.
5.  Print a message to the console indicating whether the test passed or failed.

For example, the test for `tridiagonal_matrix.m` looks like this:

```matlab
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
```
