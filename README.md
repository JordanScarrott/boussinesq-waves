# Boussinesq Wave Simulation

This repository contains a MATLAB simulation for solving the Boussinesq equations for surface water waves. The model is based on the work of Wei & Kirby (1995) on potential flow. It simulates the evolution of waves in a 2D domain with various floor profiles and initial conditions.

## Getting Started

### Prerequisites

*   MATLAB

### Installation

1.  Clone the repository:
    ```bash
    git clone https://github.com/your-username/boussinesq-waves.git
    ```
2.  Open MATLAB and navigate to the cloned repository's directory.

## Running the Simulation

To run the simulation, you need to create an instance of the `Boussinesq` class with the desired parameters, then call the `solve` and `displayMeshes` methods.

```matlab
% Example:
setup_params = [200, 0.01, 0.045, 0.45, 0.05, 0.05, 0.05, 10, 10, FloorProfile.FLAT, InitialCondition.EXPONENTIAL];
A = Boussinesq(setup_params);
A = A.solve();
A = A.displayMeshes();
```

## Parameters

The `Boussinesq` class constructor takes an array of 11 parameters:

1.  `iterations` (integer): Number of simulation iterations.
2.  `tol` (float): Tolerance for the corrector step of the numerical solver.
3.  `A0` (float): Amplitude of the driven waves.
4.  `h0` (float): Resting height of the water surface.
5.  `dx` (float): Spatial step in the x-direction.
6.  `dy` (float): Spatial step in the y-direction.
7.  `dt` (float): Time step.
8.  `real_x` (float): Real length of the domain in the x-direction.
9.  `real_y` (float): Real length of the domain in the y-direction.
10. `FloorProfile` (enum): The profile of the sea floor. See [Floor Profiles](#floor-profiles) for options.
11. `InitialCondition` (enum): The initial shape of the water surface. See [Initial Conditions](#initial-conditions) for options.

## Floor Profiles

You can specify the floor profile using the `FloorProfile` enumeration.

*   `FloorProfile.FLAT` (0): A flat, horizontal floor.
*   `FloorProfile.SINGLE_BAR` (1): A floor with a single submerged bar.

## Initial Conditions

The initial condition for the water surface can be set using the `InitialCondition` enumeration.

*   `InitialCondition.EXPONENTIAL` (0): An initial wave with an exponential decay from the center.
*   `InitialCondition.GAUSSIAN` (1): A 2D Gaussian-shaped initial wave.
*   `InitialCondition.PLANE` (2): A plane wave propagating across the domain.
*   `InitialCondition.SECH` (3): An initial wave described by a hyperbolic secant function, often used to model solitary waves.

## Features

*   **Predictor-Corrector Method:** The simulation uses a high-order Adams-Bashforth predictor and Adams-Moulton corrector scheme for time integration, providing accurate results.
*   **Boundary Conditions:** The simulation implements reflective and wavemaker boundary conditions. See `boundary_cond.m` for details.
*   **Filtering:** A 2D filter can be applied periodically to the solution to remove numerical noise. This is controlled by the `filtering` and `filter_period` properties in `Boussinesq.m`.
*   **Visualization:** The `displayMeshes` method provides an animated visualization of the wave propagation.
*   **Save/Load Parameters:** You can save the simulation parameters to an Excel file using the `saveParamData` method and load them using the static method `loadPresetFromFile`.

## Theoretical Background

The numerical model in this repository is based on the paper:

*   Wei, G., Kirby, J. T., Grilli, S. T., & Subramanya, R. (1995). A fully nonlinear Boussinesq model for surface waves. Part 1. Highly nonlinear unsteady waves. *Journal of Fluid Mechanics*, *294*, 71-92.

The paper presents a model based on fully nonlinear Boussinesq equations. This approach provides a more accurate representation of wave dynamics, especially for highly nonlinear and unsteady waves, compared to standard Boussinesq models. The model uses a high-order predictor-corrector method for time integration, which is reflected in the use of the Adams-Bashforth and Adams-Moulton schemes in this codebase. This allows for the simulation of complex wave phenomena such as shoaling and breaking with greater fidelity.

## License

This project is licensed under the terms of the LICENSE file.
