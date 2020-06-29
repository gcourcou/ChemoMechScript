import matplotlib.pyplot as plt
import numpy as np
np.pi

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
    mitotic_cell_number = [0,0,0,0,0,0,0,0,0,0]
    mitotic_frequency = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    #skip_frame+=10
    for i in range(skip_frame_start+1,len(dict_from_file["mitotic position"])-skip_frame_end):
        for j in range(0,len(dict_from_file["mitotic position"][i])):
            if dict_from_file["mitotic position"][i][j]!=[0.,0.]:
                Lap=dict_from_file["MF position"][i]-dict_from_file["L"][i][1]
                for k in range(0,10):
                  if (dict_from_file["MF position"][i] - (k+1)*Lap/10) <= dict_from_file["mitotic position"][i][j][0]:
                    if dict_from_file["mitotic position"][i][j][0] < (dict_from_file["MF position"][i] - k*Lap/10):
                        mitotic_cell_number[k] += 1
                for l in range(0,10):
                  mitotic_frequency[l] = mitotic_cell_number[l] / dict_from_file["cell number in strip"][i][l]
    
        
    plt.figure(figsize=(10,5),dpi=400)
    plt.title("mitotic_frequency")
    plt.plot([0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95],mitotic_frequency)
    plt.xticks(np.arange(0, 1, step=0.1))
    plt.xlabel("Position (x/$L_{AP}$) ",fontsize=20)
    plt.ylabel("Mitotic Frequency",fontsize=20)
    plt.savefig("mitotic_frequency_"+str(division_type)+".png",bbox_inches='tight',dpi=400)
