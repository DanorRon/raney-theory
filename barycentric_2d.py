import numpy
from utility import *
from sympy import *
from sympy.geometry import *
import pyvista as pv

perm_list = [[0,1,2], [1,0,2], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]
e_matrix = summation_matrix(perm_list[0])
t_matrix = summation_matrix(perm_list[1])

pl = pv.Plotter()

pl.show_bounds(bounds=[-5,5,-5,5,-5,5], location='all')
triangle = pv.UnstructuredGrid([3, *list(range(3))], [pv.CellType.TRIANGLE], [[1,0,0], [0,1,0], [0,0,1]])
pl.add_mesh(triangle)
pl.camera.position = (5,5,5)
pl.camera.SetParallelProjection(True)

'''
tick = 0
matrix = np.identity(3)
for i in range(5):
    if tick == 0:
        matrix = matrix @ e_matrix
        tick = 1
    else:
        matrix = matrix @ t_matrix
        tick = 0
    matrix_normalized = normalize_columns(matrix)
    corner_vecs = np.transpose(matrix_normalized)
    triangle = pv.UnstructuredGrid([3, *list(range(3))], [pv.CellType.TRIANGLE], corner_vecs)
    pl.add_mesh(triangle, color=list(np.random.choice(range(256), size=3)))
    
    print("dist0: " + str(math.dist(corner_vecs[0], corner_vecs[2])))
    print("dist1: " + str(math.dist(corner_vecs[1], corner_vecs[2])))
'''

for perm in perm_list:
    corner_vecs = np.transpose(normalize_columns(e_matrix @ summation_matrix(perm)))
    triangle = pv.UnstructuredGrid([3, *list(range(3))], [pv.CellType.TRIANGLE], corner_vecs)
    pl.add_mesh(triangle, color=list(np.random.choice(range(256), size=3)))

pl.show()