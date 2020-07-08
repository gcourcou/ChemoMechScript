import matplotlib.pyplot as plt
import numpy as np
np.pi

#def  master(file="script_out.txt",division_type="uniform") :
#    with open(file) as f:  
#        dict_from_file = eval(f.read()) 
#    
#    #dict_from_file["mitotic position"]
#    #dict_from_file["MF position"]
#    #dict_from_file["L"]
#    # skip initial frames
#    skip_frame_start=0
#    # skip final frames
#    skip_frame_end=0
#    
#    for i in range(0,len(dict_from_file["mitotic position"])):
#        if dict_from_file["MF position"][i]==0.0 and dict_from_file["MF position"][i+1]!=0.0:
#            skip_frame_start=i
#    
#    for i in range(1,len(dict_from_file["mitotic position"])-1):
#        if dict_from_file["MF position"][-i]==0.0 and dict_from_file["MF position"][-i-1]!=0.0:
#            skip_frame_end=i
#    # If we do not put this filter, then skip_frame can be equal to the whole simulation since MF=0.0 when it 
#    # reaches the anterior boundary
##    false_skips_due_to_end_filter=30
##    for i in range(0,len(dict_from_file["mitotic position"])):
##        if dict_from_file["MF position"][i]==0.0 and  i<false_skips_due_to_end_filter :
##            skip_frame=i
#        
#        #print("skipping frames " + str(skip_frame))
#    mitotic_cell_number = []
#    mitotic_frequency = []
#    number_of_slice = len(dict_from_file["cell_number_in_strip"][0])
#    for i in range(0, number_of_slice):
#        mitotic_cell_number.append(0)
#        mitotic_frequency.append(0.0)
#
#    #skip_frame+=10
#    for i in range(skip_frame_start+1,len(dict_from_file["mitotic position"])-skip_frame_end):
#        for j in range(0,len(dict_from_file["mitotic position"][i])):
#            if dict_from_file["mitotic position"][i][j]!=[0.,0.]:
#                Lap=dict_from_file["MF position"][i]-dict_from_file["L"][i][1]
#                for k in range(0,number_of_slice):
#                  if (dict_from_file["MF position"][i] - (k+1)*Lap/number_of_slice) <= dict_from_file["mitotic position"][i][j][0]:
#                    if dict_from_file["mitotic position"][i][j][0] < (dict_from_file["MF position"][i] - k*Lap/number_of_slice):
#                        mitotic_cell_number[k] += 1
#                for l in range(0,number_of_slice):
#                  mitotic_frequency[l] = mitotic_cell_number[l] / dict_from_file["cell_number_in_strip"][i][l]
#    
#        
#        plt.figure(figsize=(10,5),dpi=400)
#        plt.title("mitotic_frequency")
#        x_interval = 1/number_of_slice
#        xlist = []
#        for i in range (0, number_of_slice):
#            xlist.append((0.5+i)*x_interval)
#        
#        plt.plot(xlist,mitotic_frequency)
#        plt.xticks(np.arange(0, 1, step=x_interval))
#        plt.xlabel("Position (x/$L_{AP}$) ",fontsize=20)
#        plt.ylabel("Mitotic Frequency",fontsize=20)
#        plt.savefig("mitotic_frequency_"+str(division_type)+".png",bbox_inches='tight',dpi=400)
#    



def  master_interval(file="script_out.txt",division_type="uniform",interval=2) :
    with open(file) as f:  
        dict_from_file = eval(f.read()) 
    
    #dict_from_file["mitotic position"]
    #dict_from_file["MF position"]
    #dict_from_file["L"]

    skip_frame_start=0
    # skip final frames
    skip_frame_end=0
    
    for i in range(0,len(dict_from_file["mitotic position"])):
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
    color_saturation_c=size
    color_min_c=0
    for w in range(0,size):
        i = w * interval + skip_frame_start + 1
        for j in range(0,len(dict_from_file["mitotic position"][i])):
            if dict_from_file["mitotic position"][i][j]!=[0.,0.]:
                Lap=dict_from_file["MF position"][i]-dict_from_file["L"][i][1]
                for k in range(0,number_of_slice):
                  if (dict_from_file["MF position"][i] - (k+1)*Lap/number_of_slice) <= dict_from_file["mitotic position"][i][j][0]:
                    if dict_from_file["mitotic position"][i][j][0] < (dict_from_file["MF position"][i] - k*Lap/number_of_slice):
                        mitotic_cell_number[k] += 1
                for l in range(0,number_of_slice):
                  mitotic_frequency[l] = mitotic_cell_number[l] / dict_from_file["cell_number_in_strip"][i][l]
        #plt.hist(mitotic_position_relative_to_MF, range=[-0.2,1.1],bins=23,alpha=0.45,fill=False,density=True)
        color = plt.cm.jet((w-color_min_c)/(color_saturation_c-color_min_c))
        x_interval = 1/number_of_slice
        xlist = []
        for i in range (0, number_of_slice):
            xlist.append((0.5+i)*x_interval)
    
        plt.plot(xlist,mitotic_frequency,color=color)
        
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=color_min_c, vmax=1*color_saturation_c))
    sm.set_array([]) 
    cbar=plt.colorbar(sm)
    cbar.set_label("Time",fontsize=20)
    plt.title(division_type)
    plt.xlabel("Position (x/$L_{AP}$) ",fontsize=20)
    plt.ylabel("Mitotic Frequency",fontsize=20)
    plt.savefig("mitotic_position_"+str(division_type)+"_interval"+".png",bbox_inches='tight',dpi=400)


def  master(file="script_out.txt",division_type="uniform") :
    with open(file) as f:  
        dict_from_file = eval(f.read()) 
    
    #dict_from_file["mitotic position"]
    #dict_from_file["MF position"]
    #dict_from_file["L"]
    # skip initial frames
    skip_frame_start=0
    # skip final frames
    skip_frame_end=0
    
    for i in range(0,len(dict_from_file["mitotic position"])):
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
    for i in range(skip_frame_start+1,len(dict_from_file["mitotic position"])-skip_frame_end):
        
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
                shape_frequency[k] = dict_from_file["cell_shape_in_strip_pa"][i][k][j]/dict_from_file["cell_number_in_strip_pa"][i][k]
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
        plt.savefig("shape_frequency_"+str(division_type)+".png",bbox_inches='tight',dpi=400)
        plt.show()
        plt.close()


    
    
master(file="script_out.txt",division_type="uniform")
master_interval(file="script_out.txt",division_type="uniform",interval=2)
