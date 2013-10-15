.. manusols-4-heateqn with fenics and scipy documentation master file, created by
   sphinx-quickstart on Sun Oct 13 20:56:56 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Manufactured solutions to the 2D time-dependent heat equation
=============================================================
using sympy, fenics and scipy
-----------------------------
    The **method of manufactored solutions** is a way to test numerical code. It considers a PDE, takes an analytical solution, and computes the source term by plugging the solution into the PDE. Then, using the source term, a numerical solution is computed and compared against the analytical. See this `SANDIA report <http://prod.sandia.gov/techlib/access-control.cgi/2000/001444.pdf>`_ for a thorough introduction.
   
    We consider the time dependent heat equation. This module provides a class that has the mass and stiffness matrix as attributes and that comes with methods that provides the source term, the analytical solution projected to the discrete space, and the error between a discrete and the analytical solution. 

    It uses `sympy` for symbolic representations of the solution and the source term. 

    It uses `FEniCs` for the spatial discretization, for the projections,  and for the computation of the errors.

    The attributes, the scipy sparse matrices `M`, `A`,  represent the mass and the stiffness matrix of the problem on the **inner** nodes, i.e. that the boundary values are resolved in the source term. Thus, the time discretization is ready for the implementation of, e.g., Runge-Kutta methods for the time integration.

    Get started with `run test_impeul_int` in a python shell to plot some error of the Implicit Euler method for various time resolutions.

Contents:

.. toctree::
   :maxdepth: 2

  Code

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

