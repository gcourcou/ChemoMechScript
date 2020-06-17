import math
import random
import pandas as pd
import numpy as np
import json


# import matplotlib.pylab as plt
import ipyvolume as ipv
from tyssue import PlanarGeometry

from tyssue.geometry.planar_geometry import PlanarGeometry as geom
from tyssue.solvers.quasistatic import QSSolver
from tyssue.dynamics.planar_vertex_model import PlanarModel as model

from tyssue.config.dynamics import quasistatic_plane_spec
from tyssue import Sheet, config

from tyssue.topology.sheet_topology import cell_division

from numpy.linalg import eig, inv, pinv


def sigmoid(x, x0, k):
    return 1 / (1 + math.exp(-(x - x0) * k))


def theta(x, x0):
    if x > x0:
        return 1
    else:
        return 0


def fitEllipse(x, y):
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    D = np.hstack((x * x, x * y, y * y, x, y, np.ones_like(x)))
    S = np.dot(D.T, D)
    C = np.zeros([6, 6])
    C[0, 2] = C[2, 0] = 2
    C[1, 1] = -1
    try:
        E, V = eig(np.dot(inv(S), C))
    except:
        print("pinv")
        E, V = eig(np.dot(pinv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:, n]
    return a

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

# depriceated
def ellipse_angle_of_rotation2(a):
    b, c, d, f, g, a = a[1] / 2, a[2], a[3] / 2, a[4] / 2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi / 2
    else:
        if a > c:
            return np.arctan(2 * b / (a - c)) / 2
        else:
            return np.pi / 2 + np.arctan(2 * b / (a - c)) / 2


def elipse_division_angle(sheet, cell_index):
    local_edge_df = sheet.edge_df.loc[sheet.edge_df["face"] == cell_index]
    x = np.array([])
    y = np.array([])
    for index, row in local_edge_df.iterrows():
        x = np.append(x, row["sx"])
        y = np.append(y, row["sy"])
    a = fitEllipse(x, y)
    phi = ellipse_angle_of_rotation(a)
    axes = ellipse_axis_length(a)
    a,b=axes
    if a>=b:
        return phi
    else:
        return phi+np.pi/2


# sigmoid
def area_dep(Vd,mu,A0,A_target,A_average,alpha):
    Delta_V=Vd*(0.015+mu*(A_target-A_average)/A0)*alpha
    #print(Delta_V)
    Delta_V=max(0,Delta_V)
    return Delta_V

def area_indep(Vd,ra,alpha):
    Delta_V=Vd*0.018*ra*alpha
    #print(Delta_V)
    Delta_V=max(0,Delta_V)
    return Delta_V

def cell_Aegerter_area(sheet,Vd=1.,mu=0.04,A0=0.5,alpha=1,long_axis_div=True) :
    mitotic_shape=[0]
    mitotic_pos=[[0.,0.]]
    for index, row in sheet.face_df.iterrows():
        if row['population_variable'] == 'A':
            if row["time_for_growth"]<Vd:
                A_average=sheet.face_df['area'].mean()
                A_target=row["area"]
                sheet.face_df.at[index, "time_for_growth"] += area_dep(Vd,mu,A0, A_target ,A_average,alpha)
            else:
                mitotic_shape+=[sheet.face_df.loc[index]['num_sides']]
                mitotic_pos+=[[sheet.face_df.loc[index]['x'],sheet.face_df.loc[index]['y']]]
                sheet.face_df.at[index, "time_for_growth"]=row["time_for_growth"]/2.
                if long_axis_div == True:
                    angle_div = elipse_division_angle(sheet, index)
                else:
                    angle_div = random.random() * np.pi
                daughter = cell_division(sheet, index, geom, angle=angle_div)
    return mitotic_shape,mitotic_pos


def cell_Aegerter_uni(sheet,Vd=1.,mu=0.04,A0=0.5,alpha=1,long_axis_div=True):
    # note that input paramters are same as above for practical utility
    mitotic_shape=[0]
    mitotic_pos=[[0.,0.]]
    for index, row in sheet.face_df.iterrows():
        if row['population_variable'] == 'A':
            if row["time_for_growth"]<Vd:
                sheet.face_df.at[index, "time_for_growth"] += area_indep(Vd,row['uniform_growth_parameter'],alpha)
            else  :
                mitotic_shape+=[sheet.face_df.loc[index]['num_sides']]
                mitotic_pos+=[[sheet.face_df.loc[index]['x'],sheet.face_df.loc[index]['y']]]
                sheet.face_df.at[index, "time_for_growth"]=row["time_for_growth"]/2.
                if long_axis_div == True:
                    angle_div = elipse_division_angle(sheet, index)
                else:
                    angle_div = random.random() * np.pi
                daughter = cell_division(sheet, index, geom, angle=angle_div)
                sheet.face_df.at[index,'uniform_growth_parameter']=0.25+1.5*np.random.random()
                sheet.face_df.at[daughter,'uniform_growth_parameter']=0.25+1.5*np.random.random()
    return mitotic_shape,mitotic_pos


def cell_GS_sigmoid(
    sheet, x0=0.5, k=15, gamma_G=0.25, gamma_S=0.1, t_mech=1, long_axis_div=True
):
    mitotic_shape=[0]
    for index, row in sheet.face_df.iterrows():
        if random.random() <= sigmoid(row["area"], x0, k) * gamma_G * t_mech:
            mitotic_shape+=[ sheet.face_df.loc[index]['num_sides'] ]
            if long_axis_div == True:
                angle_div = elipse_division_angle(sheet, index)
            else:
                angle_div = random.random() * np.pi
            daughter = cell_division(sheet, index, geom, angle=angle_div)
            sheet.face_df.at[index, "cell_cycle"] = "G"
            sheet.face_df.at[daughter, "cell_cycle"] = "G"
    return mitotic_shape

def cell_GS_step(
    sheet, x0=0.5, k=15, gamma_G=0.25, gamma_S=0.1, t_mech=1, long_axis_div=True
):
    mitotic_shape=[0]
    for index, row in sheet.face_df.iterrows():
        if random.random() <= theta(row["area"], x0) * gamma_G * t_mech:
            mitotic_shape+=[ sheet.face_df.loc[index]['num_sides'] ]
            if long_axis_div == True:
                angle_div = elipse_division_angle(sheet, index)
            else:
                angle_div = random.random() * np.pi
            daughter = cell_division(sheet, index, geom, angle=angle_div)
            sheet.face_df.at[index, "cell_cycle"] = "G"
            sheet.face_df.at[daughter, "cell_cycle"] = "G"
    return mitotic_shape


def cell_GS_uniform(
    sheet, x0=0.5, k=15, gamma_G=0.25, gamma_S=0.1, t_mech=1, long_axis_div=True
):
    mitotic_shape=[0]
    for index, row in sheet.face_df.iterrows():
        if random.random() <= gamma_G * t_mech:
            mitotic_shape+=[ sheet.face_df.loc[index]['num_sides'] ]
            if long_axis_div == True:
                angle_div = elipse_division_angle(sheet, index)
            else:
                angle_div = random.random() * np.pi
            daughter = cell_division(sheet, index, geom, angle=angle_div)
            sheet.face_df.at[index, "cell_cycle"] = "G"
            sheet.face_df.at[daughter, "cell_cycle"] = "G"
    return mitotic_shape

# depreciated linear div
def cell_GS_old(sheet, amin=0.5, amax=0.6, gamma_G=0.25, gamma_S=0.1, t_mech=1):
    for index, row in sheet.face_df.iterrows():
        if row["cell_cycle"] == "G":
            # update the probability
            if row["area"] < amin:
                pass
            elif row["area"] > amax:
                if random.random() <= (amax - amin) * gamma_G * t_mech / amin:
                    sheet.face_df.at[index, "cell_cycle"] = "S"
                    # grow the cell according to geometric distribution first, then divide it.
                    p = gamma_S * t_mech
                    sheet.face_df.at[index, "time_for_growth"] = np.random.geometric(
                        p, size=1
                    )[0]
            else:
                if random.random() <= (row["area"] - amin) * gamma_G * t_mech / amin:
                    sheet.face_df.at[index, "cell_cycle"] = "S"
                    p = gamma_S * t_mech
                    sheet.face_df.at[index, "time_for_growth"] = np.random.geometric(
                        p, size=1
                    )[0]
        elif row["cell_cycle"] == "S":
            # update the probability(in fact this part never changes)
            # update the division

            if row["time_in_growth"] < row["time_for_growth"]:
                #sheet.face_df.at[index, "prefered_area"] += 1 / (
                #    sheet.face_df.at[index, "time_for_growth"]
                #)
                sheet.face_df.at[index, "time_in_growth"] += 1
            else:
                sheet.face_df.at[index, "prefered_area"] = 1
                sheet.face_df.at[index, "time_in_growth"] = 0
                angle_div = elipse_division_angle(sheet, index)
                # angle_div = random.random()*np.pi
                daughter = cell_division(sheet, index, geom, angle=angle_div)
                sheet.face_df.at[index, "cell_cycle"] = "G"
                sheet.face_df.at[daughter, "cell_cycle"] = "G"
