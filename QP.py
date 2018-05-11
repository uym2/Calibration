# Wrapper functions; copied from https://scaron.info/blog/quadratic-programming-in-python.html

import numpy
from cvxopt import matrix, spmatrix
from cvxopt.solvers import options, qp
import quadprog

def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])
    sol = qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return numpy.array(sol['x']).reshape((P.shape[1],))

def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -numpy.vstack([A, G]).T
        qp_b = -numpy.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    #return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]

def main():
    M = numpy.array([[1., 2., 0.], [-8., 3., 2.], [0., 1., 1.]])
    P = numpy.dot(M.T, M)
    q = numpy.dot(numpy.array([3., 2., 3.]), M).reshape((3,))
    G = numpy.array([[1., 2., 1.], [2., 0., 1.], [-1., 2., -1.]])
    h = numpy.array([3., 2., -2.]).reshape((3,))

    print(cvxopt_solve_qp(P, q, G, h))
    print(quadprog_solve_qp(P, q, G, h))

if __name__=="__main__":
    main()    
