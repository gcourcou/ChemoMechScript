# to make this file #plt.show() was replaced with ##plt.show()
# use agg backend matplotlib.use('Agg')
# dont forget to use sys.argv to input folder name! attempt=sys.argv[1] instead of input
import sys
sys.path.insert(0,"/Users/georgecourcoubetis/Project/Computational/Github/tyssue_git_fork")
#sys.path.insert(0,"/scratch/courcoub/tyssue/tools/tools/tyssue")

import time

start = time.time()
import pandas as pd
import numpy as np
import json
import matplotlib

#matplotlib.use("Agg")
import matplotlib.pylab as plt
import matplotlib.animation as animation
import ipyvolume as ipv
from tyssue import Sheet, config
from tyssue import PlanarGeometry
from tyssue.geometry.planar_geometry import PlanarGeometry as geom
from tyssue.solvers.quasistatic import QSSolver
from tyssue.dynamics.planar_vertex_model import PlanarModel as model
from tyssue.draw import sheet_view, sheet_view_GC_colorbar
from tyssue.stores import load_datasets
from tyssue.topology.sheet_topology import remove_face, cell_division
from scipy.spatial import Voronoi
from tyssue.config.geometry import planar_spec
from tyssue.generation import hexa_grid2d, from_2d_voronoi
#from division_functions_witharea import cell_GS
from itertools import count
from tyssue.io import hdf5

# for histogram bin number
from scipy.stats import iqr

# utils for dirs and out
import os
import sys

# data collector
script_data = {}

# data collector called for txt outputs of data
def script_data_stamp():
    try:
        os.remove("script_out.txt")
        print("script_out.txt delted")
    except:
        print("script_out.txt doesnt exist yet")
    end = time.time()
    print("time elapsed")
    total_time = end - start
    print(total_time)
    script_data["total real time"] = total_time
    script_data["cell number"] = store_cell_number
    script_data["tissue area"] = store_tissue_area
    script_data["MF position"] = store_MF_position
    script_data["Mech Timer"] = store_mech_timer
    script_data["mitotic position"] = store_mitotic_position
    script_data["L"] = store_tissue_length
    script_data["cell_number_in_strip"] = store_cell_number_in_strip
    script_data["cell_number_in_strip_pa"] = store_cell_number_in_strip_pa
    script_data["cell_shape_in_strip_pa"] = store_cell_shape_in_strip_pa
    script_data["Posterior area"] = store_P_area
    script_data["Anterior area"] = store_A_area
    f = open("script_out.txt", "w")
    f.write(str(script_data))
    f.close()
    return 


# Read parameters from file
parameters = {}
with open("parameters.txt") as f:
    for line in f:
        (key, val) = line.split()
        # print(line.split())
        try:
            parameters[str(key)] = float(val)
        except:
            parameters[str(key)] = str(val)


proliferation_type = parameters["proliferation_type"]
#long_axis_div      = parameters["long_axis_div"]
print(proliferation_type)
#print(long_axis_div)
proliferation_time_dependent = parameters["proliferation_time_dependent"]
print(proliferation_time_dependent)
# choose yourr mode

# conversion_r = pixels/micrometer
# conversion_t = s/t_mech
# Conversion is calculated with respect to original t_mech
# so we need to account if we change t_mech
og_t_mech=0.2
#1100.596923908011
#687.4683737855494 incorrect conversion
#1181.3627550441906 objective MF speed
#_5 915.3565322296259 Lposterior fit vs t_mech
#_6 874.0061531374303 Lposterior fit vs t_mech noticed that MF speed is time dependent
if proliferation_type=="area":
    from division_functions_aegerter import cell_Aegerter_area as cell_GS
    if proliferation_time_dependent=="exponential":
        def f_alpha(PL,k0=4*(10**(-5)),delta=0.0107 ,conversion_r=0.23232956642491454,conversion_t=915.3565322296259):
            conversion_t=conversion_t*parameters["conversion_t_magnitude"]
            # changing to units of current t_mech
            conversion_t=conversion_t*( parameters["t_mech"]/og_t_mech  )
            value=k0*np.exp(-1*delta*PL*(1/conversion_r))/(0.015*(1/conversion_t))
            print(value)
            return value
    elif proliferation_time_dependent=="no":
        def f_alpha(PL,k0=4*(10**(-5)),delta=0.0107 ,conversion_r=0.23232956642491454,conversion_t=915.3565322296259):
            value=parameters["proliferation_magnitude"]*( parameters["t_mech"]/og_t_mech  )
            print(value)
            return value        
elif proliferation_type=="uniform":
    from division_functions_aegerter import cell_Aegerter_uni  as cell_GS
    if proliferation_time_dependent=="exponential":
        def f_alpha(PL,k0=4*(10**(-5)),delta=0.0107,conversion_r=0.23232956642491454,conversion_t=915.3565322296259):
            conversion_t=conversion_t*parameters["conversion_t_magnitude"]
            # changing to units of current t_mech
            conversion_t=conversion_t*( parameters["t_mech"]/og_t_mech  )
            value=k0*np.exp(-1*delta*PL*(1/conversion_r))/(0.018*(1/conversion_t))
            print(value)
            return value
    elif proliferation_time_dependent=="no":
        def f_alpha(PL,k0=4*(10**(-5)),delta=0.0107 ,conversion_r=0.23232956642491454,conversion_t=915.3565322296259):
            value=parameters["proliferation_magnitude"]*( parameters["t_mech"]/og_t_mech  )
            print(value)
            return value

# output data management
#import os
#import sys

attempt = sys.argv[1]

name = "out_" + str(attempt)
os.makedirs(name, exist_ok=True)
#moved after import
#os.chdir(name)
#################


solver = QSSolver(with_collisions=False, with_t1=True, with_t3=True)

nx = parameters["nx"]
ny = parameters["ny"] + 1
pos_noise = parameters["pos_noise"]

previously_grown_eye=False
if previously_grown_eye==True:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    dsets = hdf5.load_datasets(os.path.join(__location__,"realistic_tissue.hf5") )
    #dsets = hdf5.load_datasets("./"+load_dir+'/sheet.hf5')
    specs = config.geometry.planar_sheet()
    sheet = Sheet('periodic', dsets, specs)
else:
    sheet = Sheet.planar_sheet_2d(
        "basic2D", nx=int(nx), ny=int(ny), distx=1, disty=1, noise=pos_noise
        )
# after loading tyssue i go to res dir for outputs
os.chdir(name) 
# Ellipse!
rx = parameters["rx"]
sigma = 172.0 / 74.0
ry = rx * sigma

sheet = sheet.extract_bounding_box_GC_2dellipse(rx, ry, coords=["x", "y"])


from tyssue.topology.base_topology import merge_border_edges

merge_border_edges(sheet)
#


PlanarGeometry.update_all(sheet)
sheet.get_opposite()
sheet.sanitize()
geom.update_all(sheet)
# ## Set up the model
nondim_specs = config.dynamics.quasistatic_plane_spec()
dim_model_specs = model.dimensionalize(nondim_specs)
sheet.update_specs(dim_model_specs, reset=True)

print(
    "Number of cells: {}\n"
    "          edges: {}\n"
    "          vertices: {}\n".format(sheet.Nf, sheet.Ne, sheet.Nv)
)

# ## Minimize energy
res = solver.find_energy_min(sheet, geom, model)

# ## View the result and setup display
draw_specs = config.draw.sheet_spec()
draw_specs["vert"]["visible"] = False
# size of edges
draw_specs["edge"]["head_width"] = 0.1
coords = ["x", "y"]
sheet.face_df["col"] = np.linspace(0.0, 1.0, num=sheet.face_df.shape[0])
cmap = plt.cm.get_cmap("viridis")
color_cmap = cmap(sheet.face_df.col)
draw_specs["face"]["visible"] = True
draw_specs["face"]["color"] = color_cmap
draw_specs["face"]["alpha"] = 0.5
draw_specs["axis"]={"autoscale" : False,"x_min":-30,"x_max":80,"y_min":-60,"y_max":120 }
print(draw_specs)
fig, ax = sheet_view(sheet, coords, **draw_specs)
fig.set_size_inches(7, 6)
# plt.show()
## end view

proteins = {"0": "y_concentration"}


# Initialize morphogens in dataframe
for N_p, p_name in proteins.items():
    sheet.face_df.insert(1, p_name, 0.0)

# define boundary source terms carefully since tissue moves
sheet.face_df.insert(1, "on_boundary", False)
for index, row in sheet.edge_df.iterrows():
    if row["opposite"] == -1:
        face_id = row["face"]
        sheet.face_df.at[face_id, "on_boundary"] = True
# source_term_width=parameters['source_term_width']
source_term_width = 0.2 * (sheet.face_df["x"].max() - sheet.face_df["x"].min())
with_axis_offset = sheet.face_df["x"].max() - source_term_width
sheet.face_df.insert(1, "y_boundary_source_term", False)
# x_pos_for_source = sheet.face_df['x'].max()-source_term_width
for index, row in sheet.face_df.iterrows():
    if row["x"] > with_axis_offset and row["on_boundary"]:
        sheet.face_df.at[index, "y_boundary_source_term"] = True
        local_edge_df = sheet.edge_df.loc[sheet.edge_df["face"] == index]
        for index2, row2 in local_edge_df.iterrows():
                # check if edge points t0 the boundary, otherwise target the opposite face for some boundary activation
                if row2["opposite"]!=-1:
                    # find opposite edge
                    nearest_neighbor_edge_index=row2["opposite"]
                    # find face that edge belongs to
                    nearest_neighbor_face_index=sheet.edge_df.at[nearest_neighbor_edge_index, "face"]
                    sheet.face_df.at[nearest_neighbor_face_index, "y_boundary_source_term"] = True


# MF interactions
MF_low = parameters["MF_low"]
MF_high = parameters["MF_high"]
MF_contract = parameters["MF_contract"]
MF_relax = parameters["MF_relax"]

# Steps=t_f/t_mech evalutation fo mf and div match
t_mech = parameters["t_mech"]
t_proliferation=parameters["t_proliferation"]
t_f = parameters["t_f"]
t_plot = parameters["t_plot"]

# Initialize parameters for equations
y_dif = parameters["y_dif"]
y_dec = parameters["y_dec"]
h_auto = parameters["h_auto"]

# set division params
amin = parameters["amin"]
amax = parameters["amax"]
gamma_G = parameters["gamma_G"]
gamma_S = parameters["gamma_S"]
boundary_tension = float(parameters["boundary_tension"])


#
Dif_mag = {0: y_dif}
Decay_mag = {0: y_dec}
#
auto_prod_mag = {0: h_auto, 1: []}

# auto_prod_mag[1] is updated in mech_interaction in the loop
auto_prod_mag[0] = [
    -h_auto * (int(row["y_boundary_source_term"]) * 2 - 1)
    for index, row in sheet.face_df.iterrows()
]
# Cell div pars

# Cell grow and divide params
growth_control_by = parameters["growth_control_by"]

## Initialize other cell centric variabls not included in tyssue
cell_vars = {
    "0": "cell_cycle",
    "1": "time_in_cycle",
    "2": "population _variable",
    "3": "time_in_growth",
    "4": "time_for_growth",
}
int_cell_cycle = parameters["int_cell_cycle"]
sheet.face_df["cell_cycle"] = int_cell_cycle
sheet.face_df.insert(1, "time_in_cycle", 0)

# loading file causes error so equal is better
#sheet.face_df.insert(1, "population_variable", "P")
sheet.face_df["population_variable"]=["P" for i in range(0,sheet.face_df.shape[0])]
sheet.face_df.insert(1, "time_in_growth", 0)
#sheet.face_df.insert(1, "time_for_growth", 0)

# Initialize params for aegerter growth
sheet.face_df['uniform_growth_parameter'] = 0.25+1.5*np.random.random( sheet.face_df.shape[0])
sheet.face_df[ "time_for_growth"]=np.random.random( sheet.face_df.shape[0] )

# sheet.face_df.insert(1,'on_boundary',False)

# MOVED DOWN
## Obtain protein concentration, conver to numpy array and convert to normal list for iterator
#h0 = []
#for N_p, p_name in proteins.items():
#    h0 += sheet.face_df[p_name].to_numpy().tolist()
# opposite face in the df
#sheet.get_opposite()

import numpy as np

# use len()? for faster computation
# Ncells=len(sheet.face_df)

# This function gives the diffusion contribution to dy/dt for cell with face index i
# N_p is a integer variable that dictates which protein species each function is adressed to
# Dif_mag={0:2,1:4}
def Dif(t, y, i, N_p, mag):
    # obtain all edges with desired faceid
    local_edge_df = sheet.edge_df.loc[sheet.edge_df["face"] == i]
    # print ("I AM CALLED")
    # Here I iterate through edge df to calculate diffusion terms
    term = 0.0
    for index, row in local_edge_df.iterrows():
        # if it has no oppposite face there is no diffusion
        if row["opposite"] == -1:
            term += 0.0
            # print ( "empty bc term")
        # otherwise add diffusion propotional to length
        else:
            # opposite holds value of edge not of opposite face sheet.edge_df.at[row["opposite"],'face']
            # term+= -(y[i+N_p*len(sheet.face_df)]-y[ sheet.edge_df.at[row["opposite"],'face'] +N_p*len(sheet.face_df) ])*row["length"]*mag[N_p]
            term += (
                -(
                    y[i + N_p * len(sheet.face_df)]
                    - y[
                        sheet.edge_df.at[row["opposite"], "face"]
                        + N_p * len(sheet.face_df)
                    ]
                )
                * mag[N_p]
            )
    # print (term)
    return term


# This function gives the decay contribution to dy/dy for cell with face index i
# Decay_mag={0:0.5,1:0.25}
def Decay(t, y, i, N_p, mag):
    return -mag[N_p] * y[i + N_p * len(sheet.face_df)]


# note.txt
# hh_activation= is a dict that holds

# theta function for all but controlled in such a way to include MF source cells
def Auto_Production(t, y, i, N_p, row, mag):
    return theta(y[i + N_p * len(sheet.face_df)], mag[N_p][i])


def theta(x, x0):
    if x > x0:
        return 1
    else:
        return 0.0


from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# This function gives a list of dy/dt correspondign to the y's using the contributions of diffusion,decay,reaction etc
def cell(t, y):
    fun = []
    # each cell is defined by its face index
    # 0...Ncell is the calse but lets keep it general for now
    for N_p, p_name in proteins.items():
        for index, row in sheet.face_df.iterrows():
            fun += [
                Dif(t, y, index, int(N_p), Dif_mag)
                + Decay(t, y, index, int(N_p), Decay_mag)
                + Auto_Production(t, y, index, int(N_p), row, auto_prod_mag)
            ]
        # print ("fun")
        # print (fun)
    return fun


def animate_cells2(timer, chem_name, string):
    # String tells us how to name file after image
    plt.clf()
    # sheet.face_df['col'] = np.array([ sol.y[i][-1] for i in range(0,sheet.face_df.shape[0]) ])
    sheet.face_df["col"] = sheet.face_df[chem_name]
    cmap = plt.cm.get_cmap("viridis")
    color_cmap = cmap(sheet.face_df.col)
    draw_specs["face"]["color"] = color_cmap
    draw_specs["face"]["color_bar"] = True

    #    sheet.face_df['visible'] = True
    #    for index,row in sheet.face_df.iterrows():
    #        if (row['at_y_boundary'] == True):
    #            sheet.face_df.at[index,'visible'] = False
    # fig, ax= sheet_view(sheet, coords, **draw_specs)
    fig, ax1, ax2 = sheet_view_GC_colorbar(sheet, coords, **draw_specs)
    fig.set_size_inches(6, 6)
    fig.suptitle(chem_name + " frame " + str(timer), fontsize=14)
    # fig.set_title(chem_name+' frame '+str(timer))
    plt.savefig("image" + string + chem_name + "{0:0=2d}".format(timer) + ".png",dpi=400)
    # plt.axis('off')
    # plt.show()
    plt.close()


def animate_cells_MF(timer, chem_name, string):
    # String tells us how to name file after image
    plt.clf()
    sheet.face_df["col"] = np.array(
        [
            int(row["population_variable"] == "MF") * 1.0
            for index, row in sheet.face_df.iterrows()
        ]
    )
    # sheet.face_df['col'] = sheet.face_df[chem_name]
    cmap = plt.cm.get_cmap("viridis")
    color_cmap = cmap(sheet.face_df.col)
    draw_specs["face"]["color"] = color_cmap
    draw_specs["face"]["color_bar"] = True

    #    sheet.face_df['visible'] = True
    #    for index,row in sheet.face_df.iterrows():
    #        if (row['at_y_boundary'] == True):
    #            sheet.face_df.at[index,'visible'] = False
    # fig, ax= sheet_view(sheet, coords, **draw_specs)
    fig, ax1, ax2 = sheet_view_GC_colorbar(sheet, coords, **draw_specs)
    fig.set_size_inches(6, 6)
    fig.suptitle("MF " + str(timer), fontsize=14)
    # fig.set_title(chem_name+' frame '+str(timer))
    plt.savefig("image" + string + chem_name + "{0:0=2d}".format(timer) + ".png",dpi=400)
    # plt.axis('off')
    # plt.show()
    plt.close()


def animate_cells_stateS(timer, string):
    plt.clf()
    x_ref = []
    for index, row in sheet.face_df.iterrows():
        if row["cell_cycle"] == "S":
            x_ref.append(row["x"])
    # Freedman–Diaconis rule but twice the resolution
    if len(x_ref) > 0:
        #        bin_width=2*iqr(x_ref)/(len(x_ref)**(1./3.))
        #        bin_N=int((max(x_ref)-min(x_ref))/bin_width+1)
        val = int(sheet.face_df["x"].max() - sheet.face_df["x"].min())
        bin_N = 2 * val
        plt.hist(x_ref, bins=bin_N)
        plt.title("Distribution of cells in state M in x dir frame " + str(timer))
    plt.savefig("image" + string + "{0:0=2d}".format(timer) + ".png")
    # plt.axis('off')
    # plt.show()
    plt.close()


def plot_chem(timer, chem_name, string):
    # String tells us how to name file after image
    plt.clf()
    x_y_array = np.array(
        [
            [row["x"], row["y"], row[chem_name]]
            for index, row in sheet.face_df.iterrows()
        ]
    )
    plt.plot([i[0] for i in x_y_array], [i[2] for i in x_y_array], "bx", alpha=0.25)
    x_y_ref = []
    val = (sheet.face_df["y"].max() + sheet.face_df["y"].min()) / 2.0
    for item in x_y_array:
        if val - 1.0 < item[1] < val + 1.0:
            x_y_ref.append(item)
    #plt.plot([i[0] for i in x_y_ref], [i[2] for i in x_y_ref], "rx")
    
    bin_N=10
    average_chem=[ [0] for i in range(0,bin_N)]
    Range = (sheet.face_df["x"].max() - sheet.face_df["x"].min())
    bin_size = Range/bin_N
    minimum_x = sheet.face_df["x"].min() 
    bin_x = [ minimum_x+bin_size*i for i in range(0,bin_N)  ]
    for item in x_y_array:
        for i in range(0,bin_N):
            if minimum_x+bin_size*i < item[0] <= minimum_x+bin_size*(i+1):
                average_chem[i]+= [item[2]]
    for i in range(0,bin_N):
        average_chem[i]=np.average(average_chem[i])

    plt.plot(bin_x, average_chem, "rx") 
    
    plt.title(chem_name + " vs A-P position " + str(timer) + "centered " + str(val))
    plt.savefig("image" + string + chem_name + "{0:0=2d}".format(timer) + ".png")
    # plt.axis('off')
    # plt.show()
    plt.close()


def plot_topology(timer, sheet):
    #
    # 1. frequency of number of sides
    #sheet.face_df.shape[0]
    plt.figure()
    num_of_sides_values = sheet.face_df["num_sides"].value_counts().keys().tolist()
    num_of_sides_counts = sheet.face_df["num_sides"].value_counts().tolist()
    num_of_sides_percentages=[100*item/sheet.face_df.shape[0] for item in num_of_sides_counts]
    plt.cla()
    plt.bar(num_of_sides_values, num_of_sides_percentages)
    fig.set_size_inches(12, 5)
    plt.ylabel("Pn")
    plt.xlabel("number of sides")
    # plt.axis([0,40,0,0.8])
    plt.title("cell percentages of different number of sides")
    plt.savefig("image" + "cellnumber_topology" + "{0:0=2d}".format(timer) + ".png")
    # plt.show()
    plt.close()
    # 2. average area vs number of sides
    # moved this outsdide function since it should by called once
    # sheet.face_df.insert(1,'rlarea',0.0)
    for index, row in sheet.face_df.iterrows():
        sheet.face_df.at[index, "rlarea"] = (
            sheet.face_df.at[index, "area"] / sheet.face_df["area"].mean()
        )
    rlarea_mean_dif_sides = sheet.face_df.groupby(["num_sides"])["rlarea"].mean()

    plt.cla()

    # print(rlarea_mean_dif_sides)
    plt.plot(rlarea_mean_dif_sides, marker="o")
    fig.set_size_inches(12, 5)

    plt.ylabel("<An>/<A>")
    plt.xlabel("number of sides")
    # plt.axis([0,40,0,0.8])
    plt.title("relative area of the cell versus number of sides")
    plt.savefig("image" + "cellarea_topology" + "{0:0=2d}".format(timer) + ".png")
    # plt.show()
    plt.close()


def mf_inf_cell(
    sheet, lower_limit, upper_limit, hh_inf_area, relax_area, string="y_concentration"
):
    # this function describes the cells that are constricted as a consequence of hedgehog
    # lower limit and upper limit restrict the density of hh which have influence on the cells
    # hh_inf_area describe the area of the cell that influenced by hh
    # this function is not complete, since we need to have hh_density in real situation
    # in other words, we need a relation between hh_density and x axis
    # initialize the hh density for all cells, since we do not have the density in real situation
    # this line just changes the hh_density for the forth cell, and will be changed after we get the density in real situation
    for index, row in sheet.face_df.iterrows():
        if lower_limit < row[string] < upper_limit:
            sheet.face_df.at[index, "prefered_area"] = hh_inf_area
            sheet.face_df.at[index, "cell_cycle"] = "P"
        elif row[string] >= upper_limit:
            # restore area value if mf has passed
            sheet.face_df.at[index, "prefered_area"] = relax_area
            sheet.face_df.at[index, "cell_cycle"] = "P"
    # For a certain cell, if the hh_density is between the restriction we set, the area of the cell will change to hh_inf_area
    # Find energy minimum

    # with local
    # res = solver.find_energy_min(sheet, geom, model)
    # or global
    # (this takes 4 times more time to compute)

    # res = solver.find_energy_min_GC_all(sheet, geom, model)

    print("Successfull gradient descent? ", res["success"])
    sheet.sanitize()


def visualization(i):
    animate_cells2(i, "y_concentration", "mech")
    animate_cells_MF(i, "y_concentration", "MF")
    plot_chem(i, "y_concentration", "MFvsx")
    plot_chem(i, "area", "areavsx")
    #animate_cells_stateS(i, "Sfreqvsx")
    plot_topology(i, sheet)


def mechanical_reaction(tyssue):
    mf_inf_cell(
        tyssue, MF_low, MF_high, MF_contract, MF_relax, string="y_concentration"
    )
    # define boundary source terms carefully since tissue moves
    sheet.face_df["on_boundary"] = False
    for index, row in sheet.edge_df.iterrows():
        if row["opposite"] == -1:
            face_id = row["face"]
            sheet.face_df.at[face_id, "on_boundary"] = True
    source_term_width = 0.2 * (sheet.face_df["x"].max() - sheet.face_df["x"].min())
    with_axis_offset = sheet.face_df["x"].max() - source_term_width
    sheet.face_df["y_boundary_source_term"] = False
    for index, row in sheet.face_df.iterrows():
        if row["x"] > with_axis_offset and row["on_boundary"]:
            sheet.face_df.at[index, "y_boundary_source_term"] = True
            local_edge_df = sheet.edge_df.loc[sheet.edge_df["face"] == index]
            for index2, row2 in local_edge_df.iterrows():
                # check if edge points t0 the boundary, otherwise target the opposite face for some boundary activation
                if row2["opposite"]!=-1:
                    # find opposite edge
                    nearest_neighbor_edge_index=row2["opposite"]
                    # find face that edge belongs to
                    nearest_neighbor_face_index=sheet.edge_df.at[nearest_neighbor_edge_index, "face"]
                    sheet.face_df.at[nearest_neighbor_face_index, "y_boundary_source_term"] = True
    auto_prod_mag[0] = [
        -h_auto * (int(row["y_boundary_source_term"]) * 2 - 1)
        for index, row in sheet.face_df.iterrows()
    ]
    # auto_prod_mag[1]=[0.0 for i in range(0,sheet.face_df.shape[0])]
    for index, row in sheet.face_df.iterrows():
        # auto_prod_mag[1][index]=1000
        if row["y_concentration"] < MF_low:
            sheet.face_df.at[index, "population_variable"] = "A"
        elif MF_low < row["y_concentration"] < MF_high:
            sheet.face_df.at[index, "population_variable"] = "MF"
            # set source on
            # auto_prod_mag[1][index]=z_auto
        else:
            sheet.face_df.at[index, "population_variable"] = "P"
    # z_auto

#store_mitotic_position_temp=[]
# cell_growth_and_division_2D(sheet_2D_test, alpha_c, t_M_mean, t_I_mean, m, string='c_dpp'):
def cell_grow_and_divide(tyssue):
    global store_mitotic_position
    # cell_growth_and_division_2D(tyssue, 0.5, 3, 4, 3, string='y_concentration')
    # cell_growth_and_division_2D(tyssue, 0.5, 3*2, 4*2, 3, string='z_concentration')
    #cell_GS(tyssue, amin, amax, gamma_G, gamma_S, t_mech)
    if store_MF_position[-1]==0:
        #PL=store_tissue_length[-1][0]-store_tissue_length[-1][1]
        PL=0.0
    else :
        #PL=store_MF_position[-1]-store_tissue_length[-1][1]
        PL=store_tissue_length[-1][0]-store_MF_position[-1]
    mitotic_shape, mitotic_position = cell_GS(sheet,1.,0.04,0.5,f_alpha(PL),long_axis_div=False)
    # store division position relative to MF
    store_mitotic_position+=[mitotic_position]
    tri_faces = sheet.face_df[sheet.face_df["num_sides"] < 4].index

    cells = sheet.face_df.shape[0]

    while len(tri_faces):

        remove_face(sheet, tri_faces[0])

        tri_faces = sheet.face_df[sheet.face_df["num_sides"] < 4].index
    sheet.get_opposite()
    for index, row in sheet.edge_df.iterrows():
        face_id = row["face"]
        if row["opposite"] == -1:
            sheet.face_df.at[face_id, "on_boundary"] = True
            sheet.edge_df.at[index, "line_tension"] = boundary_tension
        else:
            sheet.edge_df.at[index, "line_tension"] = 0.12
            sheet.face_df.at[face_id, "on_boundary"] = False


store_cell_number = []
store_tissue_area=[]
store_MF_position = []
store_mech_timer = []
# in cell grow and divide
store_mitotic_position = []
store_tissue_length = []
store_P_area = []
store_A_area = []

#mitotic frequency related parameters
store_cell_number_in_strip = []
store_cell_number_in_strip_pa = []
store_cell_shape_in_strip_pa = []

def data_collection(i, tyssue, cell_number, cell_number_in_strip, cell_number_in_strip_pa, cell_shape_in_strip_pa, tissue_area, MF_position, mech_timer, tissue_length,P_area, A_area):
    cell_number += [tyssue.face_df.shape[0]]
    tissue_area +=[tyssue.face_df['area'].sum()]
#    s_arr=np.array([ int(row['cell_cycle']=='S')*1.  for index, row in tyssue.face_df.iterrows() ])
#    mitosis_index +=[np.sum(s_arr) /tyssue.face_df.shape[0] ]
    x_y_array = np.array(
        [
            [row["x"], row["y"], row["population_variable"] == "MF"]
            for index, row in tyssue.face_df.iterrows()
        ]
    )
    x_y_ref = []
    val = (tyssue.face_df["y"].max() + tyssue.face_df["y"].min()) / 2.0
    for item in x_y_array:
        if val - 1.0 < item[1] < val + 1.0 and item[2] == True:
            x_y_ref.append(item[0])
    MF_mean_xpos = np.mean(x_y_ref)
    if np.isnan(MF_mean_xpos):
        MF_position += [0.0]
        MF_position_now = 0.0
    else:
        # MF position is defined as the x position of the MF cells at the midpoint of the eye disc
        MF_position += [MF_mean_xpos]
        MF_position_now = MF_mean_xpos
    mech_timer += [i * t_mech]
    #  newer way to get L 

    x_y_array = np.array(
        [
            [row["x"], row["y"], row["on_boundary"] == True]
            for index, row in tyssue.face_df.iterrows()
        ]
    )
    x_y_ref = []
    #val = (tyssue.face_df["y"].max() + tyssue.face_df["y"].min()) / 2.0
    for item in x_y_array:
        if val - 1.0 < item[1] < val + 1.0 and item[2] == True:
            x_y_ref.append(item[0])
    Lmax=max(x_y_ref)
    print("LMAX " + str(Lmax) )
    Lmin=min(x_y_ref)
    print("LMIN " + str(Lmin) )
    tissue_length += [ [Lmax,Lmin]  ]
    
    P_area_sum = 0.0
    A_area_sum = 0.0
    for index, row in sheet.face_df.iterrows():
        if row["population_variable"] == "P":
            P_area_sum += row["area"]
        elif row["population_variable"] == "A":
            A_area_sum += row["area"]
    P_area += [P_area_sum]
    A_area += [A_area_sum]
    
    cell_number_in_x_strip = []
    for k in range(0, int(parameters["number_of_slice"]) ):
        cell_number_in_x_strip.append(0)
    Lap=MF_position_now-Lmin
    strip_width = Lap/parameters["number_of_slice"]
    for index, row in sheet.face_df.iterrows():
        for j in range(0, int(parameters["number_of_slice"]) ):
            if (MF_position_now - (j+1)*strip_width) <= row['x'] and row['x'] < (MF_position_now - j*strip_width):
                 cell_number_in_x_strip[j] += 1
    cell_number_in_strip += [cell_number_in_x_strip]
    
    
    # start shape per strip analysis
    cell_shape_in_strip_pa_frame=[]
    
    cell_number_in_x_strip_pa = []
    # number of polygon classes is 6 from 4 to 10
#    init_polygon_class = list(np.zeros(6) )
#    
    for k in range(0, int(parameters["number_of_slice_pa"]) ):
        cell_number_in_x_strip_pa.append(0)
        # i want it to be [ [number_6,number_5,...], [ strip 2], [ strip 3]... ]
        cell_shape_in_strip_pa_frame.append([0,0,0,0,0,0])
    L=Lmax-Lmin
    strip_width_pa = L/parameters["number_of_slice_pa"]
    for index, row in sheet.face_df.iterrows():
        for j in range(0, int(parameters["number_of_slice_pa"]) ):
            if (Lmax - (j+1)*strip_width_pa) <= row['x'] and row['x'] < (Lmax - j*strip_width_pa):
                 cell_number_in_x_strip_pa[j] += 1
                 polygon_class=row['num_sides']
                 # minus four since zero corresponds to a four sided cell
                 try:
                     cell_shape_in_strip_pa_frame[j][polygon_class-4]+=1
                 except:
                     pass
    
    cell_number_in_strip_pa += [cell_number_in_x_strip_pa ]
    cell_shape_in_strip_pa  += [cell_shape_in_strip_pa_frame]
    # depreciated
    #tissue_length += [ [tyssue.face_df['x'].max(),tyssue.face_df['x'].min()]  ]

def chemo_mech_iterator(
    sheet,
    initial_concentration,
    derivative_function,
    mechanical_reaction,
    cell_grow_and_divide,
    visualization,
    **kwargs
):
    # mechanical_reaction should take only the tyssue object as argument and output
    # derivative function returns n-value vector args (t,x)
    # plot_function should take only the step number as argument i.e. def visualization(step) : animate_cells2(step,"mech")
    steps = int(kwargs["t_f"]/kwargs["t_mech"])
    # only proliferation steps
    steps_proliferation=int(kwargs["t_proliferation"]/kwargs["t_mech"])
    for i in range(0+steps_proliferation, steps+steps_proliferation):
        print(i)
        # time evolve by one t_mech step
        mechanical_reaction(sheet)
        # reaction diffusion
        #print(auto_prod_mag)
        sol = solve_ivp(
            derivative_function,
            [kwargs["t_mech"] * i, kwargs["t_mech"] * (i + 1)],
            initial_concentration,
        )
        sol_last_concentration = [item[-1] for item in sol.y]
        # put chem result in tissue object
        for N_p, p_name in proteins.items():
            sheet.face_df[p_name] = sol_last_concentration[
                len(sheet.face_df) * int(N_p) : len(sheet.face_df) * (int(N_p) + 1)
            ]

        # update face_df so we know where is the edge cells so we can slow them down
        for index, row in sheet.edge_df.iterrows():
            if row["opposite"] == -1:
                face_id = row["face"]
                sheet.face_df.at[face_id, "on_boundary"] = True
        
        
        data_collection(i, sheet, store_cell_number, store_cell_number_in_strip, store_cell_number_in_strip_pa, store_cell_shape_in_strip_pa,  store_tissue_area, store_MF_position, store_mech_timer, store_tissue_length, store_P_area, store_A_area)
        cell_grow_and_divide(sheet)
        solver.find_energy_min(sheet, geom, model)

        # print(sheet.face_df['time_for_M'])
        # print(sheet.face_df['prefered_area'])
        initial_concentration = []
        # necessary if mech part changes chem part
        for N_p, p_name in proteins.items():
            initial_concentration += sheet.face_df[p_name].to_numpy().tolist()
        if kwargs["plot"] == True and i%t_plot==0:
            visualization(int(i/t_plot))
        # data collection moved before cell_grow since i use the store arrays for my growth rate
        #data_collection(i, sheet, store_cell_number, store_tissue_area, store_mitosis_index, store_MF_position, store_mech_timer, store_tissue_length)
        if (i%int(kwargs["stamp"]))==0:
            script_data_stamp()
    return sol

# only cell_grow and divide
def proliferation(sheet,**kwargs):
    steps=int(kwargs["t_proliferation"]/kwargs["t_mech"])
    for i in range(0,steps):
        mechanical_reaction(sheet)
        # moved data_collection up since i use info from data for cell growthh
        data_collection(i, sheet, store_cell_number, store_cell_number_in_strip, store_cell_number_in_strip_pa, store_cell_shape_in_strip_pa, store_tissue_area, store_MF_position, store_mech_timer, store_tissue_length,store_P_area, store_A_area)
        cell_grow_and_divide(sheet)
        solver.find_energy_min(sheet, geom, model)
        if kwargs["plot"] == True and i%t_plot==0:
            visualization( int(i/t_plot) )
        # data collection moved before cell_grow since i use the store arrays for my growth rate  
        #data_collection(i, sheet, store_cell_number, store_tissue_area, store_mitosis_index, store_MF_position, store_mech_timer, store_tissue_length) 
    return

proliferation(sheet,t_proliferation=t_proliferation,t_mech=t_mech,plot=True)

# Obtain protein concentration, conver to numpy array and convert to normal list for iterator
h0 = []
for N_p, p_name in proteins.items():
    h0 += sheet.face_df[p_name].to_numpy().tolist()
sheet.get_opposite()

sol = chemo_mech_iterator(
    sheet,
    h0,
    cell,
    mechanical_reaction,
    cell_grow_and_divide,
    visualization,
    plot=True,
    stamp=20,
    t_mech=t_mech,
    t_f=t_f,
    t_proliferation=t_proliferation,
    N_protein_types=0,
)

########3 END OF SoLUTiON SCRIPT ##########
script_data_stamp()

# move back up for the sake of organization
os.chdir("..")
