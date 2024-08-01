import math
import numpy as np
import pyvista as pv
from numpy import random
import itertools
from utility import *

def many_tops(pl):
    points = []
    
    perm_list = [[0,1], [1,0]]
    base = len(perm_list)
    block_length = 11 #Only generates periodic sequences with period 13 (no less), meaning 2^13 combs
    for tick in range(base**block_length):
        periodic_block = []
        
        base_tick = to_base(tick, 2) #type str
        while len(base_tick) < block_length: #pad on the left
            base_tick = "0" + base_tick
        for digit in base_tick:
            periodic_block.append(perm_list[int(digit)])

        limcone_vertices = np.transpose(normalize_columns(generate_limcone(periodic_block, repeating=True))) #vertex 3 is the corner pt
        #presumably the two nontrivial vertices will be equal but are they never not equal?
        #Are all three vertices ever the corner pt for a nonconstant periodic perm seq?
        central_point = limcone_vertices[0]
        if is_central(central_point):
            top = minimal_top(central_point.tolist())
            normalized_top = normalize_vec(top)
            #pl.add_mesh(pv.Point((0,0,0), top), line_width=5)
            points.append(normalized_top)
        else:
            print("not central, tick = " + bin(tick)) #if not central, don't plot
    pl.add_mesh(pv.PolyData(points), color='blue')

    #print(periodic_block)
    #print(top)

def many_limcones(pl):
    points = []

    #perm_list = [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]
    #perm_list = [[0,2,1], [2,0,1]]
    perm_list = [[0,1,2], [1,0,2]] #Only plot limcones with a corner vertex
    base = len(perm_list)
    block_length = 11
    for tick in range(base**block_length):
        periodic_block = []

        base_tick = to_base(tick, 2) #type str
        while len(base_tick) < block_length: #pad on the left
            base_tick = "0" + base_tick
        for digit in base_tick:
            periodic_block.append(perm_list[int(digit)])

        #This is just any arbitrary limcone--it'd be good to restrict to only the rank 2 limcones or the r2 limcones with a corner vertex, to compare to the topped pts
        limcone_vertices = np.transpose(normalize_columns(generate_limcone(periodic_block, repeating=True)))
        face_vertex = limcone_vertices[0]
        print(str(tick) + "  " + str(face_vertex))
        points.append(face_vertex)
    pl.add_mesh(pv.PolyData(points), color='red')

pl = initialize_plotter(shape=(1,2))
pl.subplot(0,0)
many_tops(pl)
pl.subplot(0,1)
many_limcones(pl)
pl.show()

#Find angles and how they correspond to the perm seqs
#Find mapping to convert from 3d plane to 2d