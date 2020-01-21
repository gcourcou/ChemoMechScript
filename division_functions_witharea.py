import math
import random
import pandas as pd
import numpy as np
import json


#import matplotlib.pylab as plt
import ipyvolume as ipv
from tyssue import PlanarGeometry

from tyssue.geometry.planar_geometry import PlanarGeometry as geom
from tyssue.solvers.quasistatic import QSSolver
from tyssue.dynamics.planar_vertex_model import PlanarModel as model

from tyssue.config.dynamics import quasistatic_plane_spec
from tyssue import Sheet,config

from tyssue.topology.sheet_topology import cell_division

from numpy.linalg import eig, inv, pinv



def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    try:
        E, V =  eig(np.dot(inv(S), C))
    except:
        print ('pinv')
        E, V = eig(np.dot(pinv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a
def ellipse_angle_of_rotation2( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else: 
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2

def elipse_division_angle(sheet,cell_index):
    local_edge_df=sheet.edge_df.loc[sheet.edge_df['face'] == cell_index]
    x=np.array([])
    y=np.array([])
    for index, row in local_edge_df.iterrows():
        x=np.append(x,row['sx'])
        y=np.append(y,row['sy'])
    a = fitEllipse(x,y)
    phi = ellipse_angle_of_rotation2(a)
    return phi+np.pi/2.

def cell_GS(sheet,amin=0.5,amax=0.6,gamma_G=0.25,gamma_S=0.1,t_mech=1):
    for index,row in sheet.face_df.iterrows():
        if row['population_variable'] == 'A':
            if row['cell_cycle'] == 'G':
                #update the probability                                                                                                               
                if row['area'] < amin:
                    sheet.face_df.at[index,'probability_GtoS'] = 0.0
                elif row['area'] > amax:
                    #sheet.face_df.at[index,'probability_GtoS'] = (amax-amin)*gamma_G*dt/amin                                                         
                    if  random.random()<=(amax-amin)*gamma_G*t_mech/amin  :
                        #sheet.face_df.at[index,'prefered_area']=1.5                                                                                  
                        sheet.face_df.at[index,'cell_cycle'] = 'S'
                        #grow the cell according to geometric distribution first, then divide it.
                        p = gamma_S*t_mech
                        sheet.face_df.at[index,'time_for_growth'] = np.random.geometric(p,size=1)[0]
                else:
                    #print((amax-amin)*gamma_G*dt/amin)                                                                                               
                    #sheet.face_df.at[index,'probability_GtoS'] = (sheet.face_df.at[index,'area']-amin)*gamma_G*dt/amin                               
                    if  random.random() <= (row['area']-amin)*gamma_G*t_mech/amin:
                        #sheet.face_df.at[index,'prefered_area']=1.5                                                                                  
                        sheet.face_df.at[index,'cell_cycle'] = 'S'
            elif row['cell_cycle'] == 'S':
                #update the probability(in fact this part never changes)                                                                              
                #sheet.face_df.at[index,'probability_div'] = gamma_S*dt                                                                               
                #update the division
                
                if row['time_in_growth'] < row['time_for_growth'] :
                    sheet.face_df.at[index,'prefered_area'] += 1/(sheet.face_df.at[index,'time_for_growth'])
                    sheet.face_df.at[index,'time_in_growth'] += 1
                else:
                    sheet.face_df.at[index,'prefered_area'] = 1
                    sheet.face_df.at[index,'time_in_growth'] = 0
                    angle_div = elipse_division_angle(sheet,index)
                    #angle_div = random.random()*np.pi                                                                                                
                    daughter = cell_division(sheet, index, geom, angle=angle_div)
                    sheet.face_df.at[index, 'cell_cycle'] = 'G'
                    sheet.face_df.at[daughter, 'cell_cycle'] = 'G'
