import math
import numpy as np
import pyvista as pv
from numpy import random
import itertools
from utility import *
import matplotlib.pyplot as plt

def central_fan():
    points = []
    
    perm_list = [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]
    base = len(perm_list)
    block_length = 2 #Only generates periodic sequences with period 13 (no less), meaning 2^13 combs
    for tick in range(base**block_length):
        periodic_block = []
        
        base_tick = to_base(tick, base) #type str
        while len(base_tick) < block_length: #pad on the left
            base_tick = "0" + base_tick
        for digit in base_tick:
            periodic_block.append(perm_list[int(digit)])

        limcone_vertices = np.transpose(normalize_columns(generate_limcone(periodic_block, repeating=True))) #vertex 4 is the corner pt xxxxxx
        #presumably the nontrivial vertices will be equal but are they never not equal?
        #Are all the vertices ever the corner pt for a nonconstant periodic perm seq?   
        central_point = limcone_vertices[0]
        if is_central(central_point) and np.array(periodic_block)[:, 2].ptp() != 0: #don't want to count topping vectors/faces, might not be central (doubly toppable)
            top = minimal_top(central_point.tolist())
            normalized_top = normalize_vec(top)
            #pl.add_mesh(pv.Point((0,0,0), top), line_width=5)
            points.append(normalized_top)
        else:
            print("not central, tick = " + base_tick + ", sum = " + str(central_sum(central_point)), ", point = " + str(central_point)) #if not central, don't plot
    
    return points

def next_depth_fans(vertices, curr_fans, perms=list(itertools.permutations([0, 1, 2, 3]))):
    """
    Given a set of fans with a certain depth m, returns the corresponding fans with depth m + 1
    The set of fans is represented as two lists:
        vertices contains the rational vectors
        current_fans contains the fans of each rational vector at corresponding indices
    """

    next_vertices = []
    next_fans = []

    #perms_limited = [[0,1,2], [1,0,2]]
    summations = [summation_matrix(perm) for perm in perms]

    #Iterate through each rational vertex
    for i in range(len(curr_fans)):
        curr_vertex = vertices[i]
        curr_fan = curr_fans[i]
        for subdivision_matrix in summations:
            next_vertex = normalize_vec(subdivision_matrix @ curr_vertex)
            next_points = [normalize_vec(subdivision_matrix @ element) for element in curr_fan]            
            
            next_vertices.append(next_vertex)
            next_fans.append(next_points)
    
    return next_vertices, next_fans


def reduce_fan_dim(fan):
    return [barycentric_to_cartesian_3d(point) for point in fan]



    #print(periodic_block)
    #print(top)

#perms = [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]
perms = list(itertools.permutations([0, 1, 2, 3]))
perms_first = [[3,0,1,2], [3,0,2,1], [3,1,0,2], [3,1,2,0], [3,2,0,1], [3,2,1,0],
               [0,3,1,2], [0,3,2,1], [1,3,0,2], [1,3,2,0], [2,3,0,1], [2,3,1,0],
               [0,1,3,2], [0,2,3,1], [1,0,3,2], [1,2,3,0], [2,0,3,1], [2,1,3,0]]
perms_half = [[0,1,2,3], [2,3,0,1], [3,0,1,2], [3,1,2,0]] #to avoid duplicates
perms_half = perms

pl = initialize_plotter_4D(shape=(1,1))
pl.subplot(0,0)

all_vertices = []
all_fans = []

corner_point = [np.array([0,0,0,1])] #One corner
#corners = [permute(point, perm) for perm in perms for point in corner_point] #List of points; all corners
corners = [np.array([0,0,0,1]), np.array([0,0,1,0]), np.array([0,1,0,0]), np.array([1,0,0,0])]

central_fan_points = central_fan() #list of points, all at one corner

#I think I actually need the inverse permutation, see prop 2.3

central_fans = [[inverse_permute(point, perm) for point in central_fan_points] for perm in perms_half] #list of lists of points, each sublist is a fan; all corners
depth0_color = 'blue'
for fan in central_fans:
    #print(fan)
    fan_3d = reduce_fan_dim(fan)
    pl.add_mesh(pv.PolyData(fan_3d), color=depth0_color)


all_vertices = corners
all_fans = central_fans

#first subdivision fans at (0,0,1)
first_vertices_en, first_fans_en = next_depth_fans(corner_point, [central_fan_points], perms=perms_first)

#all first subdivision fans
first_vertices = [permute(point, perm) for perm in perms_half for point in first_vertices_en]
first_fans = [[permute(point, perm) for point in fan] for perm in perms_half for fan in first_fans_en] #for each fan, create a permuted fan for each perm

all_vertices.extend(first_vertices)
all_fans.extend(first_fans)

depth1_color = 'red'
for vertex in first_vertices:
    vertex_3d = barycentric_to_cartesian_3d(vertex)
    pl.add_mesh(pv.PolyData(vertex_3d), color='black')
for fan in first_fans:
    fan_3d = reduce_fan_dim(fan)
    pl.add_mesh(pv.PolyData(fan_3d), color=depth1_color)


#central_fans = [central_fan_points]

next_vertices = first_vertices
next_fans = first_fans
for i in range(1):
    next_vertices, next_fans = next_depth_fans(next_vertices, next_fans)

    all_vertices.extend(next_vertices)
    all_fans.extend(next_fans)

    depth_color = list(np.random.choice(range(256), size=3))

    for vertex in next_vertices:
        vertex_3d = barycentric_to_cartesian_3d(vertex)
        pl.add_mesh(pv.PolyData(vertex_3d), color='black')
    

    for fan in next_fans:
        fan_3d = reduce_fan_dim(fan)
        pl.add_mesh(pv.PolyData(fan_3d), color=depth_color)

pl.show()

#Find angles and how they correspond to the perm seqs
#Find mapping to convert from 3d plane to 2d

'''
e = summation_matrix([2,0,1])
b = central_fan_points[0]
print(b)
print(normalize_vec(e @ b))
'''

#where is each corner vertex mapped to by each matrix?
#How do I even know that a matrix multiplied by a limcone is a limcone?
#Multiplying M by a vector doesn't always increase the depth

#etetet...: the point
#eetetet: new point
#What kind of limcone do I get? The vertex seems to have depth 0 so does the interior also have depth 0? I think so
#The issue is, this should be a depth 0 limcone (I think) but it's generated with the depth 1 limcones in this process
#And I can't really just remove it, it's a limcone that otherwise wouldn't be counted.
#Is there a theorem about the depth of an ultimately convergent sequence? Like if it's periodic after a certain point what's the depth

#tetetet: the other central fan point. I should not have the propagation include this because it has already been counted in depth 0

#Any other starting value should be a limcone with depth 1

#Need to account for this at the other corners too

#So I don't think stating the perm as eetetet... gives the whole story. The issue is rather that the depth of a perm is the value k at which P_k(n) is thereafter constant.
#Thus, if the perm that I stack onto the beginning has the same P(3) (which for e_3 is 3) as the rest of the sequence, the depth doesn't increase.
#So I should exclude those two permutations from the perms which I multiply by to get the next depth.
#I should draw out all the central fan mappings to make sure every corner of the first subdivision has something mapped to it still (even with this exclusion).
#Does this rule apply for subsequent subdivisions? Will any perms not increase the depth if I use them then?
#I think it just applied for the first subdivision: Consider a summation matrix multiplied by a depth 1 vector; must increase the depth.

#print(permute([0,1,2], [2,1,0]))

#I think the vanes that have been found here are the limcones accepting perm seqs with ultimately constant P_k(n). As the first m perms are given for a certain depth m, 
#in the limit as m --> infinity, I essentially get every combination (maybe)

#Ofc I'm only considering the periodic sequences here, so this also doesn't include the non-periodic limcones (with constant P_k(n) or not)