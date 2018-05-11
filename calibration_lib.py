from dendropy import Tree
import numpy as np
from QP import quadprog_solve_qp, cvxopt_solve_qp
from math import exp,log, sqrt
from scipy.stats.stats import pearsonr
from scipy.optimize import minimize

def deviation_from_clock(tree):
    factor = []
    for node in tree.postorder_node_iter():
        if node is tree.seed_node:
            continue
        f = np.random.gamma(exp(1.5),scale=1/exp(1.5))
        #f = np.random.lognormal(1,0.4)
        factor.append(f)
        node.edge_length = node.edge_length*f

    return factor    


def calibrate_log_opt(tree,smpl_times):
    def f(x):
        #print(x)
        return sum([log(y)*log(y) for y in x[:-1]])
        #return sum([1.0/y*1.0/y-2*1.0/y for y in x[:-1]])

    def g(x,a):    
        return sum([ x[i]*a[i] for i in range(len(x)) ])


    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    cons = []
    
    idx = 0
    
    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1
        if node.is_leaf():
            node.constraint = [0.0]*(N+1)
            node.constraint[node.idx] = node.edge_length
            node.constraint[N] = -smpl_times[node.taxon.label]
        else:
            children = list(node.child_node_iter())
                       
            a = np.array([ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N+1) ])
            cons.append({'type':'eq','fun':g,'args':(a,)})

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length

    x0 = [1.]*N + [0.006]
    bounds = [(0.00000001,999999)]*(N+1)

    result = minimize(f,x0,bounds=bounds,constraints=cons,method="SLSQP")
    x = result.x
    #f = f0 
    s = x[N]
    
    print("Clock rate: " + str(s))
    
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            #node.time = node.smplTime
            node.time = -node.constraint[N]
            #print(node.taxon.label + " " + str(node.time))
        else:
            children = list(node.child_node_iter())
            node.time = children[0].time - x[children[0].idx]*children[0].edge_length/s
            #print(node.label +  " " + str(node.time))
        #if node is not tree.seed_node:
        #    node.edge_length *= x[node.idx]/s
    
    return f

def calibrate_with_sampling_time(tree,smpl_times,verbose=False):
    if verbose:
        print("Start calibrating tree. Setting up the matrices and constraints...")

    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2


    h = np.array([0.0]*(N+1)).reshape((N+1,))
    b = np.array([0.0]*(n-1)).reshape((n-1,)) 
    #b = np.array([0.0]*(n-1)).reshape((n-1,)) 
    #b = []

    G = np.negative(np.identity(N+1))
    #P = N*np.identity(N) - np.ones((N,N))
    #q = np.array([0.0]*N).reshape((N,))

    P = np.identity(N+1)
    P[N][N] = 0
    q = np.array([-2.0]*N+[0]).reshape((N+1,)) 
    
    #A = [1.0]*N
    A = np.zeros((n-1,N+1))
    #A[0] = [1.0]*N + [0]
    #A = None

    idx = 0
    r = 0
    c = 0
    L = 1000.0

    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1

         
        if node is not tree.seed_node: 
            P[node.idx][node.idx] = node.edge_length*node.edge_length/(node.edge_length+c/L)
            q[node.idx] = -2.0*node.edge_length*node.edge_length/(node.edge_length+c/L)
        

        if node.is_leaf():
            node.constraint = [0.0]*(N+1)
            node.constraint[node.idx] = node.edge_length
            node.constraint[N] = -smpl_times[node.taxon.label]
            #node.smplTime = smpl_times[node.taxon.label]
        else:
            children = list(node.child_node_iter())
                       
            a = [ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N+1) ]

            #A = np.vstack([A,a])
            #A = np.vstack([A,a]) if A is not None else a
            #b.append(children[0].smplTime - children[1].smplTime)
            A[r] = a
            r += 1

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length
                #node.smplTime = children[0].smplTime

    if verbose:
        print("Done with setting up. Started quadratic optimization ...")

    f = cvxopt_solve_qp(P,q,G,h,A,b)
    #f = quadprog_solve_qp(P,q,G,h,A,b)
    w = f[N]
    print("Clock rate: " + str(w))
    if verbose:
        print("Optimal solution found. Calibrating tree ...")
    
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            #node.time = node.smplTime
            node.time = -node.constraint[N]
            #print(node.taxon.label + " " + str(node.time))
        else:
            children = list(node.child_node_iter())
            node.time = children[0].time - f[children[0].idx]*children[0].edge_length/w
            #print(node.label +  " " + str(node.time))
        #if node is not tree.seed_node:
        #    node.edge_length *= f[node.idx]/s
    
    return f


def calibrate_tree(tree,verbose=False):
    if verbose:
        print("Start calibrating tree. Setting up the matrices and constraints...")

    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2

    h = np.array([0.0]*N).reshape((N,))
    #b = np.array([N]+[0.0]*(n-1)).reshape((n,)) 
    b = np.array([0.0]*(n-1)).reshape((n-1,)) 

    G = np.negative(np.identity(N))
    P = np.identity(N)
    #q = np.array([0.0]*N).reshape((N,))
    q = np.array([-2.0]*N).reshape((N,)) 
    
    #A = [1.0]*N
    #A = np.zeros((n,N))
    #A[0] = [1.0]*N
    A = None

    idx = 0
    r = 1
    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1
        if node.is_leaf():
            node.constraint = [0.0]*N
            node.constraint[node.idx] = node.edge_length
        else:
            children = list(node.child_node_iter())
                       
            a = [ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N) ]
            #A = np.vstack([A,a])
            A = np.vstack([A,a]) if A is not None else a
            #A[r] = a
            r += 1
            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length

    if verbose:
        print("Done with setting up. Started quadratic optimization ...")

    #f = cvxopt_solve_qp(P,q,G,h,A,b)
    f = quadprog_solve_qp(P,q,G,h,A,b)
    if verbose:
        print("Optimal solution found. Calibrating tree ...")

    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            node.edge_length *= f[node.idx]

    return f

def main():
    from sys import argv

    tree = Tree.get_from_path(argv[1],"newick")
    
    '''
    smpl_times = {}
    with open(argv[2],"r") as fin:
        fin.readline()
        for line in fin:
            name,time = line.split()
            smpl_times[name] = float(time)
    '''
    f = deviation_from_clock(tree)

    tree.write_to_path(argv[2],"newick")
   
    '''
    m=sum([1/x for x in f])/len(f)
    with open('f.txt','w') as fout:
        for x in f:
            fout.write(str(1/x/m) + "\n")
    '''    

    #f = calibrate_with_sampling_time(tree,smpl_times)
    #f = calibrate_tree(tree)
    #print(f)
    '''
    m = sum(f1)/len(f1)
    with open('f1.txt','w') as fout:
        for x in f1:
            fout.write(str(x/m) + "\n")
    '''

    #tree.write_to_path(argv[2],"newick")
    
    #print(np.corrcoef([1/x for x in f],f1))

if __name__=="__main__":
    main()     
