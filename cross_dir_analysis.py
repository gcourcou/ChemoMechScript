import os
import matplotlib
matplotlib.use("Agg")

import matplotlib.pylab as plt     
plt.rcParams['figure.figsize'] = (7,5)
font = {'family' : 'normal',
        'size'   : 20}
plt.rc('font', **font)
calc=2760/5
plt.rcParams['figure.dpi']=calc
plt.rcParams['savefig.bbox']='tight'
import numpy as np

from import_analyze_script_out import analyze


def heatmap(a,x_data,y_data):
    plt.imshow(a, cmap='hot', interpolation='nearest',origin='lower')

    #x-ticks
    nx = x_data.shape[0]
    no_labels = nx # how many labels to see on axis x
    step_x = int(nx / (no_labels - 1)) # step between consecutive labels
    x_positions = np.arange(0,nx,step_x) # pixel count at label position
    x_labels = x_data[::step_x] # labels you want to see
    x_labels=np.around(x_labels,decimals=1)
    plt.xticks(x_positions, x_labels)
    plt.xlabel('a0')
    
    #y-ticks
    
    
    ny = y_data.shape[0]
    no_labels = ny # how many labels to see on axis x
    step_y = int(ny / (no_labels - 1)) # step between consecutive labels
    y_positions = np.arange(0,ny,step_y) # pixel count at label position
    y_labels = y_data[::step_y] # labels you want to see
    y_labels=np.around(y_labels,decimals=1)
    plt.yticks(y_positions, y_labels)
    plt.ylabel('k')
    
    plt.colorbar()
    #cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    plt.savefig("phase_diagram.png")
    #plt.show()

#only difference is that bottom iso out_gen2 rather than default setting in analyze

targets=[]
top_dir=os.getcwd()

CUTOFF_DEPTH=1
for root, dirs, files in os.walk(".", topdown=False):
   if root.count(os.sep) >= CUTOFF_DEPTH:
       del dirs[:]
   for name in files:
      1==1
      #print(os.path.join(root, name))
   for name in dirs:
      print("dirs:")
      print(dirs)
      if root.count(os.sep) == CUTOFF_DEPTH-1 and name != ".git":
          #print(os.path.join(root, name))
          targets+=[os.path.join(root, name)[1:]]
print("targets:")
print(targets)
store={}
parameter_in_question="Repeats"
for directory in targets:
    if directory[-1]!="_":
        os.chdir(top_dir+directory)
        print(str(directory))
        temp=analyze(bottom="out_0/")
        # touple
        param=np.around(temp['parameters'][parameter_in_question],decimals=2)  
        if parameter_in_question == "MF_contract":
            store[np.around(1-param, decimals=2)] = temp
        else:
            store[param]=temp
        #print(os.listdir())

colorbar_name = {}
parameters = {}
with open("parameters.txt") as f:
    for line in f:
        (key, val) = line.split()
        try:
            parameters[str(key)] = float(val)
        except:
            parameters[str(key)] = str(val)
    
for key in parameters:
    if key == "MF_contract":
        colorbar_name[key] = "Apical constriction factor"
    else:
        colorbar_name[key] = key

#print(store)
print("store keys:")
print(store.keys())

#return to top dir
os.chdir(top_dir)

#import numpy as np
#import matplotlib.pyplot as plt
#xlist = np.linspace(0, 5, 5)
#ylist = np.linspace(0, 5, 5)
#X, Y = np.meshgrid(xlist, ylist)
#Z = np.sqrt(X**2 + Y**2)
#fig,ax=plt.subplots(1,1)
#cp = ax.contourf(X, Y, Z)
#fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
##ax.set_xlabel('x (cm)')
#ax.set_ylabel('A0 ')
#ax.set_xlabel('k ')
#plt.show()

# zero was not valid
#x_data = np.arange(.1,.9,.1)
x_data=[0.0 + 0.1*(i) for i in range(0,7)]  
#x_data = np.arange(1.12,1.32,.02)
x_data=np.around(x_data,decimals=2)

# colour_data defined as whether development finished or not
colour_data={}
for item_x in x_data:
    if store[item_x]['finished']:
        colour_data[item_x]='blue'
    else:
        colour_data[item_x]='red'
# mag_colour_data defined as magnitude of variable which is changed 
mag_colour_data={}
color_saturation_c=max(x_data)
color_min_c=min(x_data)
for item_x in x_data:
    mag_colour_data[item_x] = plt.cm.jet((item_x-color_min_c)/(color_saturation_c-color_min_c) )

#y_data = np.arange(5,35,5)
#y_data=np.around(y_data,decimals=1)

dict_y_x_labels={
'MF speed':[colorbar_name[parameter_in_question],"MF speed $(μm h^{-1})$"],'MF linearity':[colorbar_name[parameter_in_question],"MF linearity $R^2$"],
'final average area':[colorbar_name[parameter_in_question],"Final average area ($μm^2$)"],'final area':[colorbar_name[parameter_in_question],"Final area($μm^2$)"],'total_mitotic_number':[colorbar_name[parameter_in_question], "Total mitotic number"], 
'average_anterior_cell_area':[colorbar_name[parameter_in_question],"Average anterior cell area ($μm^2$)"],'average_posterior_cell_area':[colorbar_name[parameter_in_question],"Average posterior cell area ($μm^2$)"],
'Lp':["Time (h)","Lp ($μm$)"],'La':["Time (h)","La ($μm$)"],
'area':["Time (h)","Area ($μm^2$)"], 'anterior_area':["Time (h)","Anterior area ($μm^2$)"], 'posterior_area':["Time (h)","Posterior area ($μm^2$)"],
'average_number_of_sides_in_MF':["Time (h)","Average number of sides in MF"], 'average_area_in_MF':["Time (h)","Average area in MF ($μm^2$)"], 'MF_shape':["Time (h)","MF shape (root square deviation)"],
'anterior_cell_area':["Time (h)","Anterior cell area ($μm^2$)"],'posterior_cell_area':["Time (h)","Posterior cell area ($μm^2$)"],
    'cell number':["Time (h)","Cell number"],'cell_death':["Time (h)","Apoptotic Cell Number"],'cell_division':["Time (h)","Mitotic Cell Number"],
'mitotic_number_list':["Normalized Time","Division Number"],
'mitotic_frequency':["relative position to MF", "frequency"], 'final cell number':[colorbar_name[parameter_in_question],"Cell Number"]
}


plot_keys=['MF speed','MF linearity','final average area','final area','average_anterior_cell_area','average_posterior_cell_area','total_mitotic_number','final cell number']

a={}
for key in plot_keys:
    a[key]=[]


for key in plot_keys:
    temp=[]
    for item_x in x_data:
        temp+=[store[item_x][key]]    
    a[key]+=[temp][0]

## plot a time depnendent quanity for many regions in one plot
plot_keys_vectors=['Lp','La','area','anterior_area','posterior_area','anterior_cell_area','posterior_cell_area','cell number','cell_division',
                   'cell_death']

# plot for a range that correspodns to average t_mechh
# i use median

plot_keys_MF=['average_number_of_sides_in_MF', 'average_area_in_MF', 'MF_shape', 'mitotic_frequency']

# valid for plot_key_vectors Chi you need to add the rest for plot_key_MF and use them for nice plot
#dict_y_x_labels={'Lp':["Time (h)","($μm$)"],'La':["Time (h)","($μm$)"],
#'area':["Time (h)","($μm^2$)"], 'anterior_area':["Time (h)","($μm^2$)"], 'posterior_area':["Time (h)","($μm^2$)"]
#}
for key in plot_keys_vectors:
    plt.figure()
    # store length of array with time dep quantity
    time_dep_len=[]
    max_y=[]
    min_y=[]
    max_x=[]
    min_x=[]
    for item_x in x_data:
        #915 is the converion_t
        # accounting for t_mech
        og_t_mech=0.2
        time_array=np.arange(0,len(store[item_x][key]))*store[item_x]['parameters']['conversion_t_magnitude']*915.3565322296259*(store[item_x]['parameters']['t_mech']/og_t_mech)*0.000277778
        plt.plot(time_array,store[item_x][key],color=mag_colour_data[item_x])
        #time_dep_len+=[len(store[item_x][key])]
        max_y+=[np.max(store[item_x][key])]
        min_y+=[np.min(store[item_x][key])]
        max_x+=[np.max(time_array)]
        min_x+=[np.min(time_array)]
    x_range_up=np.max(max_x)
    x_range_down=np.min(min_x)
    plt.xlim(x_range_down,x_range_up)
    y_range_up=np.max(max_y)
    y_range_down=np.min(min_y)
    plt.ylim(y_range_down,y_range_up)
    
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=color_min_c, vmax=1*color_saturation_c))
    cbar=plt.colorbar(sm)
    cbar.set_label(colorbar_name[parameter_in_question],fontsize=20)
#    plt.title(str(key))
    print(str(key) + " time dep plot x range is " )
    print(x_range_up)
    print(x_range_down)
    # y and x labes dict_y_x_labels holds a list of labels [xlabel,ylabel]
    plt.xlabel(dict_y_x_labels[key][0])
    plt.ylabel(dict_y_x_labels[key][1])
    plt.savefig("cross_plot_"+str(key)+".png")
    plt.close()
## end time dep plot

for key in plot_keys_MF:
    plt.figure()
    # store length of array with time dep quantity
    time_dep_len=[]
    max_y=[]
    min_y=[]
    max_x=[]
    min_x=[]
    for item_x in x_data:
        #915 is the converion_t
        # accounting for t_mech
        og_t_mech=0.2
        time_array=np.arange(0,len(store[item_x][key]))*store[item_x]['parameters']['conversion_t_magnitude']*915.3565322296259*(store[item_x]['parameters']['t_mech']/og_t_mech)*0.000277778
        plt.plot(time_array,store[item_x][key],color=mag_colour_data[item_x])
        #time_dep_len+=[len(store[item_x][key])]
        max_y+=[np.max(store[item_x][key])]
        min_y+=[np.min(store[item_x][key])]
        max_x+=[np.max(time_array)]
        min_x+=[np.min(time_array)]
    x_range_up=np.max(max_x)
    x_range_down=np.min(min_x)
    plt.xlim(x_range_down,x_range_up)
    y_range_up=np.max(max_y)
    y_range_down=np.min(min_y)
    plt.ylim(y_range_down,y_range_up)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=color_min_c, vmax=1*color_saturation_c))
    cbar=plt.colorbar(sm)
    cbar.set_label(colorbar_name[parameter_in_question],fontsize=20)
#    plt.title(str(key))
    print(str(key) + " time dep plot x range is " )
    print(x_range_up)
    print(x_range_down)
    # x y labbesl from dict
    plt.xlabel(dict_y_x_labels[key][0])
    plt.ylabel(dict_y_x_labels[key][1])
    plt.savefig("cross_plot_"+str(key)+".png")
    plt.close()

for key in plot_keys:
    plt.figure()
    print(x_data)
    print(a[key])
    print(key)
    print("average" + str(np.average(a[key])) )
    print("variance" + str(np.var(a[key]) ) )
    for item_x in x_data:
        plt.scatter(item_x,store[item_x][key],color=colour_data[item_x])
    #plt.plot(x_data,a[key],'.')
#    plt.title(str(key))
    # x y labels from dict
    plt.xlabel(dict_y_x_labels[key][0])
    plt.ylabel(dict_y_x_labels[key][1])
    plt.savefig("cross_plot_"+str(key)+".png")
    plt.close()

plt.figure()
max_y = []
min_y = []

key = "mitotic_number_list"
time_key = "normalized_time_list"
for item_x in x_data:
    plt.plot(store[item_x][time_key],store[item_x][key],color=mag_colour_data[item_x])
    max_y+=[np.max(store[item_x][key])]
    min_y+=[np.min(store[item_x][key])]
x_range_up=1.0
x_range_down=0.0
plt.xlim(x_range_down,x_range_up)
y_range_up=np.max(max_y)
y_range_down=np.min(min_y)
plt.ylim(y_range_down,y_range_up)
sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=color_min_c, vmax=1*color_saturation_c))
cbar=plt.colorbar(sm)
cbar.set_label(colorbar_name[parameter_in_question],fontsize=20)
#plt.title(str(key))
# x y labbesl from dict
plt.xlabel(dict_y_x_labels[key][0])
plt.ylabel(dict_y_x_labels[key][1])
plt.savefig("cross_plot_"+str(key)+".png")
plt.close()
#test = np.random.random((6, 5))
#
#heatmap(test,x_data,y_data)
