#!/usr/bin/env python
# coding: utf-8

# In[386]:

import itertools
from sage.plot.plot3d.shapes import Shape
from sage.plot.plot3d.shapes import LineSegment
from sage.plot.plot3d.shapes import Sphere
from sage.plot.plot3d.shapes import polygons3d

rubix_dict = {'r': [1,0,0], 'l':[-1,0,0],
              'f': [0,1,0], 'b':[0,-1,0],
              'u': [0,0,1], 'd':[0,0,-1]}

def rubix_to_face(rubix_list):
    coords = []
    for full in rubix_list:
        a,b,c = 0,0,0
        for single in full:
            a += rubix_dict[single][0]
            b += rubix_dict[single][1]
            c += rubix_dict[single][2]
        coords.append([a,b,c])
    return coords


def face_to_cartesian(rubix_list):
    face_list = []
    for face in rubix_list:
        squarelist = [[i] if i != 0 else [-1,1] for i in face]
        new_face = list(itertools.product(*squarelist))
        if len(new_face) == 2:
            new_face += [tuple(2*i for i in j) for j in new_face]
        new_face = [new_face[1], new_face[0] , new_face[2], new_face[3]]
        face_list.append(new_face)
    return face_list


def palette(index):
    color_list = ['red','orange','yellow','green','lightblue','blue','darkblue','purple']
    if index < len(color_list):
        return color_list[index]
    else:
        return black

    
def show_it(rubix_coords,rotation=[0,0]):
    color_index = 0
    z = -1 * rotation[0]
    x = -1 * rotation[1]  
    face_coords = rubix_to_face(rubix_coords)
    cartesian_coords = face_to_cartesian(face_coords)
    P =  Sphere(.1, color='red')
    P += Sphere(.1, color='blue').translate(0,.01,0)
    for k in cartesian_coords:
        P += polygons3d([range(4)], k, color='gray',opacity=0.7)
#         P += polygons3d([range(4)], k, color='gray',opacity=0.7).rotate((0,0,1),z*math.pi/2).rotate((1,0,0),x*math.pi/2).scale(0.5,0.5,0.5)
        if k[0][0]/k[3][0] == k[0][1]/k[3][1]:
            P += LineSegment(k[2], k[3], 2, color=palette(color_index))
            P += LineSegment(k[2], k[3], 4, color=palette(color_index)).rotate((0,0,1),z*math.pi/2).rotate((1,0,0),x*math.pi/2).scale(0.5,0.5,0.5)
            color_index +=1
    P.show()

# O PATH : 3/3 COMPLETE
show_it(['uf','ur','ub','ul'])
show_it(['uf','ur','ub','ul','f','r','l','d'],[2,1])
show_it(['uf','ur','ub','ul','d','u'],[0,2])


# L PATH : 4/5 COMPLETE
show_it(['ul','fl','dl','db','rb','ub'],[0,0])
show_it(['ul','fl','dl','db','rb','ub','f','b'],[-1,0])
# show_it(['ul','fl','dl','db','rb','ub'],[0,3]) # FAIL
show_it(['ul','fl','dl','db','rb','ub','d','f'],[1,-1])
show_it(['ul','fl','dl','db','rb','ub','u','d'],[2,0])


# U PATH : 2/4 COMPLETE
show_it(['ul','fl','dl','db','dr','fr','ur','ub'],[0,0])
# show_it(['ul','fl','dl','db','dr','fr','ur','ub'],[1,0]) # FAIL
show_it(['ul','fl','dl','db','dr','fr','ur','ub','f','b'],[2,0])
# show_it(['ul','fl','dl','db','dr','fr','ur','ub'],[-1,0]) # FAIL

# OO PATH : 2/2 COMPLETE
show_it(['uf','ur','ub','ul','df','dr','db','dl'])
show_it(['uf','ur','ub','ul','df','dr','db','dl','r','l'],[0,1])
get_ipython().run_line_magic('pinfo', 'more')

# EQ PATH : 2/2 COMPLETE
show_it(['uf','ul','bl','db','dr','rf'])
show_it(['uf','ul','bl','db','dr','rf','u','d'],[2,0])
# A petrie hexagon

