.. manusols-4-heateqn with fenics and scipy documentation master file, created by
   sphinx-quickstart on Sun Oct 13 20:56:56 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Manufactored solutions to the 2D time-dependent heat equation
=============================================================
using sympy, fenics and scipy
-----------------------------
    This module can be used to test time integration schemes, that are implemented, e.g., in scipy, via the technique of manufactured solutions. This means that one considers a PDE (here the 2D time dependent heat equation), makes up a solution, computes a spatial discretization and the associated right hand side, and uses this right hand side to numerically compute solutions to the actual system. Thus, the solution is known and one can test the scheme for convergence and performance.

    We use `sympy` to provide symbolic representations of the solution and the right hand sides. 

    We use `FEniCs` for the spatial discretization and computation of the errors.

    And this module, that provides the coefficient matrices in `scipy`'s csr matrix format and right hand sides ready for use in time integration schemes.

    An instance of the main class as attributes the scipy sparse matrices `M`, `A` representing the mass and the stiffness matrix of the problem on the **inner** nodes and methods that return the right hand side of the discrete system, the exact solution, and the error at a given time instance (and for a given approximate solution).

Contents:

.. toctree::
   :maxdepth: 2

  Code

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

