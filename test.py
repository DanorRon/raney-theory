import math
import numpy as np
import pyvista as pv
from numpy import random
import itertools
from utility import *

perm_list = [[0,1,2], [1,0,2], [2,0,1]]
e_matrix = summation_matrix(perm_list[0])
t_matrix = summation_matrix(perm_list[1])
a_matrix = summation_matrix(perm_list[2])

'''
matrix = np.identity(3)
for i in range(18):
    rand = np.random.random_sample()
    if rand > 0.5:
        matrix = matrix @ e_matrix
    else:
        matrix = matrix @ t_matrix
    print(normalize_columns(matrix))

#for i in range(1000):
    #matrix = matrix @ a_matrix
    #print(normalize_columns(matrix))

matrix = normalize_columns(matrix)

print(matrix)
print(is_central(np.transpose(matrix)[0]))
'''
'''
print(summation_matrix([1,2,0]))
print(summation_matrix([2,1,0]))
a = normalize_columns(generate_limcone([[0,1,2], [0,2,1]], repeating=True))
print(minimal_top(np.transpose(a)[0].tolist()))
'''
'''
perm_block = [[0,1,2,3], [1,0,2,3], [0,1,3,2], [1,0,3,2]]
topped_perm_block = [[0,1,2,3,4], [1,0,2,3,4], [0,1,3,2,4], [1,0,3,2,4]]
print(generate_limcone(perm_block))
limcone_4d = normalize_columns(generate_limcone(perm_block, repeating=True))
print(limcone_4d)
print(generate_limcone(topped_perm_block))
limcone_5d = normalize_columns(generate_limcone(topped_perm_block, repeating=True))
print(limcone_5d)
print(normalize_vec(minimal_top(np.transpose(limcone_4d)[0].tolist())))
print(normalize_vec(minimal_top(np.transpose(limcone_4d)[2].tolist())))
on_line = 0.7*np.transpose(limcone_4d)[0] + 0.3*np.transpose(limcone_4d)[2]
print(normalize_vec(minimal_top(on_line.tolist())))
'''

#print(summation_matrix([0,1,2]))

#print(generate_limcone([[0,1,2], [1,0,2]], repeating=True))

'''
tick = 0
matrix = np.identity(3)
for i in range(100):
    if tick == 0:
        matrix = matrix @ e_matrix
        tick = 1
    else:
        matrix = matrix @ t_matrix
        tick = 0
    matrix_normalized = normalize_columns(matrix)
    corner_vecs = np.transpose(matrix_normalized)
    print("dist0: " + str(math.dist(corner_vecs[0], corner_vecs[2])))
    print("dist1: " + str(math.dist(corner_vecs[1], corner_vecs[2])))
'''

#print(np.transpose(normalize_columns(generate_limcone(perms=[[0,1,2], [1,2,0]], repeating=True))))

a = find_perm_seq([1,1,0])
print(a)
print(generate_limcone([[1,2,0]]))


#I think nonperiodic sequences can still create central limcones, the requirement is that the norm is ultimately >= 1. Look into this
#This might mean that everything is either central or a topping