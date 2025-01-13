import numpy as np
import scipy as sp
from scipy.optimize import minimize


def updateAllElements(fes, req, dx):
    for fe in fes:
        fe.updateState(req, dx)


def solveWithSciPyMinimize(energyFun, x0, jacFun=None, hessFun=None, callBackFun=None):
    if jacFun == None:
        return minimize(
            energyFun,
            x0=x0,
            options={"disp": True, "gtol": 1e-3},
            tol=1e-14,
            callback=callBackFun,
        )
    elif hessFun == None:
        return minimize(
            energyFun,
            x0=x0,
            jac=jacFun,
            options={"disp": True, "gtol": 1e-3},
            tol=1e-14,
            callback=callBackFun,
        )
    else:
        return minimize(
            energyFun,
            method="trust-constr",
            x0=x0,
            jac=jacFun,
            hess=hessFun,
            options={"disp": True},
            callback=callBackFun,
        )


# Using derivative free version root, as the others do not accept a callback
def solveWithSciPyRoot(fun, x0, callBackFun=None):
    if callBackFun != None:
        return sp.optimize.root(
            fun, x0=x0, tol=1e-10, method="krylov", callback=callBackFun
        )


def solveWithNewtonRaphson(gradFun, hessFun, assembler, fes, req):
    maxiter = 100
    abs_tolerance = 1e-14
    d = np.zeros(assembler.reducedSize())
    for k in range(maxiter):
        R = gradFun(d)
        K = hessFun(d)
        r_norm = np.linalg.norm(R)

        deltad = sp.linalg.solve(K, R)

        updateAllElements(fes, req, assembler.createFullVector(deltad))
        d -= deltad

        if r_norm < abs_tolerance or k > maxiter:
            break

    success = r_norm < abs_tolerance
    return success, d
