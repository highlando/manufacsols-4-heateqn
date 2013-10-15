from dolfin import *
import sympy as smp
import numpy as np
import scipy.sparse as sps
import krypy

from matplotlib import pyplot as pp

parameters.linear_algebra_backend = "uBLAS"

def tran_heat_eq2d(omega=1):
    """ An example time varying 2D manufactured solution

    :param omega: Parameter for the frequency. Defaults to 1

    :return: A symbolic expression of a scalar function depending on the space variables `x` (internal name is `x[0]`) and `y` (internal name is `x[1]`) and the time variable `t`

    """
    x, y, t = smp.symbols('x[0], x[1], t')
    time_dep = smp.cos(omega*t)
    symb_sol = time_dep * smp.exp(x)*x*x*2*y*(1-y)*(2*y-1)
    return symb_sol

def test_impeul_int():
    """ Example for the test of convergence of Implicit Euler
    """

    symb_sol = tran_heat_eq2d(omega=9)
    mesh = UnitSquareMesh(10, 10)
    V = FunctionSpace(mesh, 'CG', 2)
    prob = ProbFenicsNumpiHeat2d(symb_sol, V)

    Ntslist = [16, 32, 64]
    t0, tE = 0.0, 1.0

    u_old = prob.cursolprjtd(t0)
    err_old = prob.comp_cur_error(u_old, t0)
    err_int, err_int_list = 0.0, []

    def _sol_linsys(A, rhs, x0=None, tol=None, maxiter=None):
        """Solution of the linear system is wrapped here

        for easy modification
        """

        return np.atleast_2d(sps.linalg.spsolve(A, rhs)).T
        # return krypy.linsys.gmres(A, rhs, x0=x0,
        #                             tol=tol,
        #                             maxiter=maxiter
        #                           )['xk']

    print 'Number of timesteps, global error'
    for Nts in Ntslist:
        tau = (tE-t0)/Nts
        curl_A = prob.M - tau * prob.A
        for tcur in np.linspace(t0+tau, tE, Nts-1):

            f_st = prob.M * u_old + tau * prob.currhs(tcur)
            u_new = _sol_linsys(curl_A, f_st, x0=u_old)

            err_new = prob.comp_cur_error(u_new, t0)

            # pw. trap rule for err int in time
            err_int = tau * 0.5 * (err_new + err_old)
            err_old, u_old = err_new, u_new

        err_int_list.append(err_int)
        print Nts, err_int

    pp.loglog(Ntslist, err_int_list, '*', color='blue')
    pp.xlabel('Number of timesteps')
    pp.ylabel('||u-u_h||')
    pp.show()


class ProbFenicsNumpiHeat2d(object):
    """
    This class takes a sympy expression as the manufactured solution 
    :math:`u(t;x,y)` to the nonstationary heat equation

    .. math::

        u_t - \\nabla \\cdot (\\kappa \\nabla u) = f

    completed with **Dirichlet** boundary conditions and an initial
    condition

    and a FEM function space where the equation is approximated on.

    :param symbsol: a symbolic expression with parameters `t,x,y` with internal names `x[0], x[1], t`
    :param V: FeNiCs functionspace on a 2D mesh
    :param kappa: diffusion coefficient: may be a constant value or a fenics expression in `x[0]` and `x[1]`. Defaults to 1


    """
    def __init__(self, symbsol, V, kappa=1):

        x, y, t = smp.symbols('x[0], x[1], t')
        rhs = symbsol.diff(t) \
              + (kappa * (symbsol.diff(x))).diff(x) \
              + (kappa * (symbsol.diff(y))).diff(y) 

        from sympy.printing import ccode
        self.sol = Expression(ccode(symbsol), t=0.0)
        self.rhs = Expression(ccode(rhs), t=0.0)
        self.V = V

        v = TestFunction(V)
        u = TrialFunction(V)
        
        # Assemble system
        M = assemble(inner(u, v)*dx)
        A = assemble((-1) * kappa * inner(grad(u), grad(v)) * dx)
        # Convert DOLFIN representation to numpy arrays
        rows, cols, values = M.data()
        self.M = sps.csr_matrix((values, cols, rows))
        """csr matrix for the mass"""
        rows, cols, values = A.data()
        self.A = sps.csr_matrix((values, cols, rows))
        """csr matrix representing the weak discrete 
        :math:`\\nabla \\cdot (\\kappa \\nabla )` operator"""
            

        # treatment of the boundary values
        nv = self.A.shape[0]
        auxu = np.zeros((nv,1))
        self.bcinds = []
        self.bc = DirichletBC(self.V, self.sol, 'on_boundary')
        self.bcdict = self.bc.get_boundary_values()
        auxu[self.bcdict.keys(),0] = self.bcdict.values()
        self.bcinds.extend(self.bcdict.keys())

        self.rhsbc = - self.A*auxu    
    
        # indices of the innernodes
        self.invinds = np.setdiff1d(range(nv),self.bcinds).astype(np.int32)
        # condense the coeff mats to the inner nodes 
        self.M = self.M[self.invinds,:][:,self.invinds]
        self.A = self.A[self.invinds,:][:,self.invinds]

        self.rhsbc = self.rhsbc[self.invinds,:]

        self.bcvals = auxu[self.bcinds]

    def currhs(self, tcur):
        """Compute the current source term 
        
        containing the actual source of the system plus the 
        contribution from the Dirichlet boundary data

        :param tcur: current time

        :return: Vector of the current rhs

        """
        self.rhs.t = tcur
        v = TestFunction(self.V)
        currhs = assemble(inner(self.rhs, v)*dx)
        currhs = currhs.array().reshape(len(currhs), 1)
        return currhs[self.invinds,:] + self.rhsbc

    def cursolprjtd(self, tcur):
        """ Get the current solution projected onto the 
        FEM space

        :param tcur: current time

        :return: vector of the current projected solution

        """
        self.sol.t = tcur
        cursolprj = project(self.sol, self.V, 
                solver_type = "lu", bcs=self.bc)
        cursolprj = cursolprj.vector()
        cursolprj = cursolprj.array().reshape(len(cursolprj), 1)
        return cursolprj[self.invinds,:]

    def expand_vec_to_fun(self, uvec):
        """expand uvec to the dolfin function representation

        :param uvec: vector of a current solution approximation

        :return: fenics function containing the values of `uvec` and the Dirichlet data

        """

        vfun = Function(self.V)
        ve = np.zeros((self.V.dim(),1))

        bcdict = self.bc.get_boundary_values()
        ve[bcdict.keys(),0] = bcdict.values()

        ve[self.invinds] = uvec 
        vfun.vector().set_local(ve)

        return vfun

    def comp_cur_error(self, uvec, tcur):
        """get the current error

        :param uvec: vector of solution approximation
        :param tcur: current time 

        :return: the approximation error 

        """

        vfun = self.expand_vec_to_fun(uvec)
        self.sol.t = tcur

        return errornorm(self.sol, vfun)


if __name__ == '__main__':
    test_impeul_int()

