import matplotlib.pyplot as plt
import numpy as np
from numpy import nan
np.pi

parameters={}
with open("parameters.txt") as f:
    for line in f:
        (key, val) = line.split()
        try:
            parameters[str(key)] = float(val)
        except:
            parameters[str(key)] = str(val)
amin = parameters["amin"]
amax = parameters["amax"]
gamma_G = parameters["gamma_G"]
gamma_S = parameters["gamma_S"]
boundary_tension = float(parameters["boundary_tension"])
t_mech=float(parameters["t_mech"])
og_t_mech=0.2
proliferation_magnitude = float(parameters["proliferation_magnitude"])
conversion_t_magnitude = float(parameters["conversion_t_magnitude"])


conversion_r=0.23232956642491454
#pixel/microm
conversion_t=915.3565322296259*conversion_t_magnitude*(t_mech/og_t_mech)
#s/t_mech
conversion_t_hr=conversion_t*0.000277778

bottom = "out_0/"

def  mitotic_frequency_sum(dict_from_file,division_type="uniform",interval=2) :
#    with open(bottom+file) as f:  
#        dict_from_file = eval(f.read()) 
    
    #dict_from_file["mitotic position"]
    #dict_from_file["MF position"]
    #dict_from_file["L"]

    skip_frame_start=0
    # skip final frames
    skip_frame_end=0
    
    for i in range(0,len(dict_from_file["mitotic position"])-1):
        if dict_from_file["MF position"][i]==0.0 and dict_from_file["MF position"][i+1]!=0.0:
            skip_frame_start=i
    
    for i in range(1,len(dict_from_file["mitotic position"])-1):
        if dict_from_file["MF position"][-i]==0.0 and dict_from_file["MF position"][-i-1]!=0.0:
            skip_frame_end=i
        #skip_frame+=20
        #print("skipping frames " + str(skip_frame))
    #skip_frame+=10
    
    mitotic_cell_number = []
    mitotic_frequency = []
    number_of_slice = len(dict_from_file["cell_number_in_strip"][0])
    for i in range(0, number_of_slice):
        mitotic_cell_number.append(0)
        mitotic_frequency.append(0.0)

    plt.figure(figsize=(10,5),dpi=400)
    size=int((len(dict_from_file["mitotic position"])-skip_frame_start-skip_frame_end)/interval-1)
    print(size)
    color_saturation_c=size*conversion_t_hr
    color_min_c=0
    out_frame = int(size/2)
    out_list = []
    for w in range(0,size):
        ylist = []
        for num in range(0, number_of_slice):
            ylist.append(0.0)
        for m in range(0, interval):
            i = w * interval + skip_frame_start + 1 + m
            for num in range(0, number_of_slice):
                mitotic_cell_number[num] = 0
                mitotic_frequency[num] = 0.0
            for j in range(0,len(dict_from_file["mitotic position"][i])):
                if dict_from_file["mitotic position"][i][j]!=[0.,0.]:
    #                print("mitotic number:")
    #                print(len(dict_from_file["mitotic position"][i]))
                    Lam=dict_from_file["MF position"][i]-dict_from_file["L"][i][1]
                    for k in range(0,number_of_slice):
                      if (dict_from_file["MF position"][i] - (k+1)*Lam/number_of_slice) <= dict_from_file["mitotic position"][i][j][0]:
                        if dict_from_file["mitotic position"][i][j][0] < (dict_from_file["MF position"][i] - k*Lam/number_of_slice):
                            mitotic_cell_number[k] += 1
    #        print(i)
    #        print("mitotic_cell_number is:")
    #        print(mitotic_cell_number)
            for l in range(0,number_of_slice):
                if dict_from_file["cell_number_in_strip"][i][l] != 0:
                    mitotic_frequency[l] = mitotic_cell_number[l] / dict_from_file["cell_number_in_strip"][i][l]
                else:
                    mitotic_frequency[l] = 0.0
            for num in range(0, number_of_slice):
                ylist[num] += mitotic_frequency[num]
#        print("cell number in strip is:")
#        print(dict_from_file["cell_number_in_strip"][i])
        #plt.hist(mitotic_position_relative_to_MF, range=[-0.2,1.1],bins=23,alpha=0.45,fill=False,density=True)
        color = plt.cm.jet((w*conversion_t_hr-color_min_c)/(color_saturation_c-color_min_c))
        x_interval = 1/number_of_slice
        xlist = []
        for i in range (0, number_of_slice):
            xlist.append((0.5+i)*x_interval)
    
        plt.plot(xlist,ylist,color=color)
        
        if w == out_frame:
            out_list = ylist
        
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=color_min_c, vmax=1*color_saturation_c))
    sm.set_array([]) 
    cbar=plt.colorbar(sm)
    cbar.set_label("Time",fontsize=20)
    plt.title(division_type)
    plt.xlabel("Position (x/$L_{AP}$) ",fontsize=20)
    plt.ylabel("Mitotic Frequency sum",fontsize=20)
    plt.savefig("mitotic_position_sum__"+str(division_type)+"_interval"+".png",bbox_inches='tight',dpi=400)
    return out_list


def  shape_frequency(dict_from_file,division_type="uniform", interval = 2) :
#    with open("out_0/"+file) as f:  
#        dict_from_file = eval(f.read()) 
    
    #dict_from_file["mitotic position"]
    #dict_from_file["MF position"]
    #dict_from_file["L"]
    # skip initial frames
    skip_frame_start=0
    # skip final frames
    skip_frame_end=0
    
    for i in range(0,len(dict_from_file["mitotic position"])-1):
        if dict_from_file["MF position"][i]==0.0 and dict_from_file["MF position"][i+1]!=0.0:
            skip_frame_start=i
    
    for i in range(1,len(dict_from_file["mitotic position"])-1):
        if dict_from_file["MF position"][-i]==0.0 and dict_from_file["MF position"][-i-1]!=0.0:
            skip_frame_end=i
    # If we do not put this filter, then skip_frame can be equal to the whole simulation since MF=0.0 when it 
    # reaches the anterior boundary
#    false_skips_due_to_end_filter=30
#    for i in range(0,len(dict_from_file["mitotic position"])):
#        if dict_from_file["MF position"][i]==0.0 and  i<false_skips_due_to_end_filter :
#            skip_frame=i
        
        #print("skipping frames " + str(skip_frame))
    init_polygon_class = len(dict_from_file["cell_shape_in_strip_pa"][0][0])
    shape_frequency = []
    number_of_slice_pa = len(dict_from_file["cell_number_in_strip_pa"][0])
    for i in range(0, number_of_slice_pa):
        shape_frequency.append(0.0)

    #skip_frame+=10
    size=int((len(dict_from_file["mitotic position"])-skip_frame_start-skip_frame_end)/interval-1)
    for m in range(0,size):
        i = m * interval + skip_frame_start + 1
        plt.figure(figsize=(10,5),dpi=400)
        size = init_polygon_class
        color_saturation_c=size
        color_min_c=0
        
        shape_frequency = []
        for w in range(0, number_of_slice_pa):
            shape_frequency.append(0.0)
        for j in range(0,init_polygon_class):
            polygon_class = j+4
            for k in range(0, number_of_slice_pa):
                if dict_from_file["cell_number_in_strip_pa"][i][k] != 0:
                    shape_frequency[k] = dict_from_file["cell_shape_in_strip_pa"][i][k][j]/dict_from_file["cell_number_in_strip_pa"][i][k]
                else:
                    shape_frequency[k] = 0.0
            x_interval_pa = 1/number_of_slice_pa
            xlist = []
            for l in range (0, number_of_slice_pa):
                xlist.append((0.5+l)*x_interval_pa)
        
            color = plt.cm.jet((j-color_min_c)/(color_saturation_c-color_min_c))
            plt.plot(xlist,shape_frequency,color=color)
            
        sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=color_min_c, vmax=1*color_saturation_c))
        sm.set_array([]) 
        cbar=plt.colorbar(sm)
        cbar.set_label("num_of_sides",fontsize=20)
        cbar.ax.set_yticklabels(['4','5','6','7','8','9'])
        plt.title(division_type)
        plt.xlabel("Position (x/$L_{AP}$) ",fontsize=20)
        plt.ylabel("shape Frequency",fontsize=20)
        plt.savefig("shape_frequency_"+str(division_type)+str(i)+".png",bbox_inches='tight',dpi=400)
        plt.close()

def mitotic_number(dict_from_file):
    
    skip_frame_start=0
    # skip final frames
    skip_frame_end=0
    
    for i in range(0,len(dict_from_file["mitotic position"])-1):
        if dict_from_file["MF position"][i]==0.0 and dict_from_file["MF position"][i+1]!=0.0:
            skip_frame_start=i
    
    for i in range(1,len(dict_from_file["mitotic position"])-1):
        if dict_from_file["MF position"][-i]==0.0 and dict_from_file["MF position"][-i-1]!=0.0:
            skip_frame_end=i
    
    total_frame = len(dict_from_file["mitotic position"])-skip_frame_start-skip_frame_end  
    
    normalized_time_list = []
    mitotic_number_list = []

    for i in range(0, total_frame):
        normalized_time = i/total_frame
        mitotic_number = len(dict_from_file["mitotic position"][i])
        normalized_time_list.append(normalized_time)
        mitotic_number_list.append(mitotic_number)
    
    plt.figure()
    plt.plot(normalized_time_list, mitotic_number_list)
    plt.title("mitotic_number")
    plt.xlabel("normalized time", fontsize=20)
    plt.ylabel("mitotic number", fontsize=20)
    plt.savefig("mitotic_number_vs_time.png")
    plt.close()
    
    return sum(mitotic_number_list)
    
    
    
    
#mitotic_frequency_sum(dict_from_file,division_type="uniform", interval=20)
#shape_frequency(dict_from_file,division_type="uniform",interval=20)
