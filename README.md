pde-solvers
===========

This repository contains some of the big c++ research codes I wrote in my past life as a scientist at LLNL. These versions are no longer actively maintained.

Each directory contains a discretized partial differential equation (PDE) solver for a particular problem arising in fluid dynamics, cell biology and finance. 

Many of the codes depend on the [Overture framework](http://www.overtureframework.org/) and would require some work to build with the current releases. However, the 1-d test cases have no dependencies and should build on a Unix box. The purpose of the 1-d test codes was to demonstrate and study some of the physics and core numerical issues arising in the full 3-D solver.

###Summary
- [CellWave](https://github.com/petrifast/pde-solvers/tree/master/CellWave):  A solver for a reaction-diffusion system of equations modeling the dynamics of spatial biochemical reactions in 1-d to 3-d model systems. 

  The 1-d version of CellWave is used for studying the interaction of diffusion and reaction and allows tuning model parameters. The N-dimensional version (N=2 or 3) uses [Overture](http://www.overtureframework.org/) and allows simulations in complex geometry described by the overset grid files.

