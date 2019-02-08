# axi-mushy-layer
Simple model for convection in an axisymmetric mushy-layer. 

Author: Jamie Parkinson (jamie.parkinson@gmail.com)

# About 
`axisymmetricMushyLayer.m` is the heart of the code. It sets up the solve and then calls `steadyStateSolver` to find the steady states.

The poisson equation for the stream function is solved by `poissonSolverTwoDimensionsAxisymm`, and the heat equation is solved by `heatSolverADI`.

Boundary conditions are specified in `mushyLayerTemperatureBoundaryConditions` and `mushyLayerPsiBoundaryConditions`.

The initial state is defined in `mushyLayerInitialState`.

Other codes in this folder:
* `extrapolateSolution` computes the solution in the channel
* `calculateThetaInfinity` computes the far field temperature for some solution
* `extrapolate_a` predicts the new chimney position, I think
* `calculateQ` computes the  heat advection  $\mathbf{U} \cdot \nabla T$
* `axiVelocitiesFromPsi` computes the velocity from the axisymmetric streamfunction.

Other subfolders:
* `cartesian` contains some cartesian versions of this code
* `util` contains useful scripts for finding data files with certain parameters, returning the properties of the mesh, and computing some diagnostics
* `benchmarking` contains scripts for testing our code against various benchmark problems
* `docs` contains the report describing the physics behind the problem, and some results from this code

