import math
import numpy as np
import pyvista as pv
from numpy import random
import itertools
from utility import *

def many_tops(pl):
    points = []
    
    perm_list = [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]
    base = len(perm_list)
    block_length = 5 #Only generates periodic sequences with period 13 (no less), meaning 2^13 combs
    for tick in range(base**block_length):
        periodic_block = []
        
        base_tick = to_base(tick, base) #type str
        while len(base_tick) < block_length: #pad on the left
            base_tick = "0" + base_tick
        for digit in base_tick:
            periodic_block.append(perm_list[int(digit)])

        limcone_vertices = np.transpose(normalize_columns(generate_limcone(periodic_block, repeating=True))) #vertex 4 is the corner pt
        #presumably the two nontrivial vertices will be equal but are they never not equal?
        #Are all three vertices ever the corner pt for a nonconstant periodic perm seq?
        central_point = limcone_vertices[0] #TODO this might not work for higher dims, check this (the first 3 columns might not be equal)
        if is_central(central_point):
            top = minimal_top(central_point.tolist())
            normalized_top = normalize_vec(top)
            #pl.add_mesh(pv.Point((0,0,0), top), line_width=5)
            points.append(normalized_top)
            print(normalized_top) #TODO write comments for dims of vectors
        else:
            print("not central, tick = " + bin(tick)) #if not central, don't plot
        #TODO map the 4d vectors onto a 3d tetrahedron before adding to points
    pl.add_mesh(pv.PolyData(points), color='blue')

    #print(periodic_block)
    #print(top)

def many_limcones(pl):
    points = []

    perm_list = [[0,1,2,3], [0,2,1,3], [1,0,2,3], [1,2,0,3], [2,0,1,3], [2,1,0,3]] #Only plot limcones with a corner vertex
    base = len(perm_list)
    block_length = 5
    for tick in range(base**block_length):
        periodic_block = []

        base_tick = to_base(tick, 6) #type str
        while len(base_tick) < block_length: #pad on the left
            base_tick = "0" + base_tick
        for digit in base_tick:
            periodic_block.append(perm_list[int(digit)])

        #This is just any arbitrary limcone--it'd be good to restrict to only the rank 2 limcones or the r2 limcones with a corner vertex, to compare to the topped pts
        limcone_vertices = np.transpose(normalize_columns(generate_limcone(periodic_block, repeating=True)))
        print(limcone_vertices)
        face_vertex = limcone_vertices[0]
        print(str(tick) + "  " + str(face_vertex))
        points.append(face_vertex)
        #TODO map the 4d vectors onto a 3d tetrahedron before adding to points
    pl.add_mesh(pv.PolyData(points), color='red')

pl = initialize_plotter_4D(shape=(1,2))
#pl.subplot(0,0)
#many_tops(pl)
pl.subplot(0,1)
many_limcones(pl)
pl.show()

#Find angles and how they correspond to the perm seqs
#Find mapping to convert from 3d plane to 2d