import numpy as np
import math

def permutation_matrix(perm):
    """
    Return a permutation matrix given a permutation written as the second row in 2-row notation
    """
    dim = len(perm)
    matrix = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            if j == perm[i]:
                matrix[i][j] = 1
    return matrix

def summation_matrix(perm):
    """
    Returns the summation matrix M_p given the permutation perm
    """
    dim = len(perm)
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            if i >= j:
                A[i][j] = 1
    B_p = permutation_matrix(perm)
    M_p = np.linalg.inv(B_p) @ A @ B_p
    return M_p

def isotone_perm(x):
    """
    Finds the permutation that when applied to x, makes x isotone (where applying is the operation x -> xP)
    Requires x to be nonnegative
    Equivalent to finding P_k+1 given x_k
    """
    x_copy = x.copy()
    perm = np.zeros((len(x)))
    for i in range(len(x)): #repeats len(x) times
        max = np.argmax(x_copy)
        perm[len(x) - i - 1] = max
        x_copy[max] = -1 # sets to value out of domain
        #print(str(max) + "  " + str(x) + "  " + str(perm))
    return perm

def next_xk(xk):
    """
    Returns x_k+1 given x_k
    Uses the matrix ideas in the proof of Proposition 2.8
    """
    P_k1 = isotone_perm(xk) #P_k+1
    M_P_k1 = summation_matrix(P_k1) # M(P_k+1)
    #x_k = M(P_k+1) @ x_k+1
    x_k1 = np.linalg.inv(M_P_k1) @ xk
    #print(str(M_P_k1) + "  " + str(x_k1))
    return x_k1

x_i = [1, math.sqrt(2)] #x0
print(x_i)
for i in range(2):
    x_i = next_xk(x_i)
    print(x_i)