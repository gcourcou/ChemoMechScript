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

from scipy.stats import chisquare

from sklearn.metrics import r2_score
from numpy import nan

from mitotic_frequency import mitotic_frequency_sum
from mitotic_frequency import shape_frequency
from mitotic_frequency import mitotic_number
#scipy.stats.chisquare(f_obs, f_exp=None, ddof=0, axis=0)
#### Input parameters and define functions for estimations

#def calculate_p_G():
#    return (np.average(dict_from_file["cell ave area"])-amin)*gamma_G*t_mech/amin
#
#def calculate_p_S():
#    return gamma_S*t_mech
#####

def r_2_and_poly(x_data,y_data,order=1):
        res=np.polyfit(x_data, y_data, order,full=True)
        #SStot=np.array(y_data)-np.average(y_data)
        #SStot=np.sum(SStot**2)
        
        #SSres= y_data-(np.poly1d(res[0])(x_data))
        #SSres= np.sum(SSres**2)
        #Rsquared=(SStot-SSres)/SStot
        Rsquared=r2_score(y_data, np.poly1d(res[0])(x_data) )
        return res[0], Rsquared
def  scatter_and_fit(target_x,target_y,text_x="Time (h)",text_y="Lp (μm)",text="none"):
    #fig=plt.figure(figsize=(7,5),dpi=calc)
    plt.figure()
    ax=plt.axes()
    #target_x=np.linspace(1,len(target_y), num=len(target_y))-1
    #t=np.linspace(0,len(target_y),1)
    #ax.set_title(text)
    plt.plot(target_x,target_y,'.')
    plt.ylabel(text_y)
    plt.xlabel(text_x)    
    res,Rsquared= r_2_and_poly(target_x,target_y)
    plt.plot(target_x, np.poly1d(res)(target_x))
    ax.text(0.5,0.8,"R^2 = "+str(Rsquared)[0:5],transform=ax.transAxes )
    #plt.ylabel("Displacement  (mm) ")
    #plt.xlabel("Spot Area ($mm^2$) ")
    #plt.savefig('cross_region_trajectory_distance_travelled_vs_displacement.png')
    #plt.savefig(text+".png",bbox_inches='tight',calc=2760/5)
    plt.savefig('Lp_vs_time_fit.png')
    #plt.show()
    plt.close()
    return res,Rsquared

def analyze(bottom="./"):
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
    with open(bottom+"script_out.txt") as f:  
        dict_from_file = eval(f.read()) 
    
    
    conversion_r=0.23232956642491454
    #pixel/microm
    conversion_t=915.3565322296259*conversion_t_magnitude*(t_mech/og_t_mech)
    #s/t_mech
    conversion_t_hr=conversion_t*0.000277778
    #hr/t_mech
    time_array=np.arange(0,len(dict_from_file["tissue area"]))*conversion_t_hr
    
    
    mitotic_frequency_list = mitotic_frequency_sum(dict_from_file,division_type="uniform", interval=20)
    
#    shape_frequency(dict_from_file,division_type="uniform",interval=20)
    
    normalized_time_list,mitotic_number_list,total_mitotic_number = mitotic_number(dict_from_file, interval = 3)
    
    plot_area_keys=["tissue area","Posterior area","Anterior area"]
    #names used for plots
    plot_area_names=["Total area","Posterior area","Anterior area"]
    plot_number_of_sides_in_MF=["average_number_of_sides_in_MF", "MF_shape"]
    plot_area_MF=["average_area_in_MF"]

    
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
    plot_frame = len(dict_from_file["mitotic position"])-skip_frame_end
    plot_time_array = []
    for j in range(0, plot_frame):
            plot_time_array.append(time_array[j])
    
    for i in range(0,len(plot_area_keys)):
        key=plot_area_keys[i]
        yname=plot_area_names[i]
        plt.figure()
        plt.ylabel(yname+" ($μm^2$)")
        plt.xlabel("Time (h)")
        plot_area_array = []
        area=np.array(dict_from_file[key])/((conversion_r)**2)
        for j in range(0, plot_frame):
            plot_area_array.append(area[j])
        plt.plot(plot_time_array,plot_area_array,'.')
        plt.savefig(yname+"_vs_time.png")
        plt.close()
    for i in range(0,len(plot_number_of_sides_in_MF)):
        key=plot_number_of_sides_in_MF[i]
        yname=key
        plt.figure()
        plt.ylabel(yname)
        plt.xlabel("Time (h)")
        ylist = []
        time_array_1 = []
        for j in range(0, len(time_array)):
            if dict_from_file[key][j] != 0.0 and np.isnan(dict_from_file[key][j]) != True:
                ylist += [dict_from_file[key][j]]
                time_array_1 += [time_array[j]]
        plt.plot(time_array_1,ylist,'.')
        plt.savefig(yname+"_vs_time.png")
        plt.close()
    for i in range(0,len(plot_area_MF)):
        key=plot_area_MF[i]
        yname=key
        plt.figure()
        plt.ylabel(yname+" ($μm^2$)")
        plt.xlabel("Time (h)")
        area= []
        time_array_1 = []
        for j in range(0, len(time_array)):
            if dict_from_file[key][j] != 0.0 and np.isnan(dict_from_file[key][j]) != True:
                area += [dict_from_file[key][j]/((conversion_r)**2)]
                time_array_1 += [time_array[j]]
        plt.plot(time_array_1,area,'.')
        plt.savefig(yname+"_vs_time.png")
        plt.close()
#    area=np.array(dict_from_file["tissue area"])/((conversion_r)**2)
#    plt.plot(area,'.')
#    plt.savefig('area_vs_time.png')
#    plt.close()
    plt.figure()
    dict_from_file["cell number"]
    plt.ylabel("Cell number ")
    plt.xlabel("Time (h)")
    cell_number_array = []
    for i in range(0, plot_frame):
        cell_number_array.append(dict_from_file["cell number"][i])
    plt.plot(plot_time_array,cell_number_array,'.')
    #plt.show()
    plt.savefig('cell_number_vs_time.png')
    plt.close()
    
    average_area=np.array(dict_from_file["tissue area"])/np.array(dict_from_file["cell number"])
    average_area=average_area/((conversion_r)**2)
    plt.figure()
    plt.ylabel("Average cell area ($μm^2$)")
    plt.xlabel("Time (h)")
    average_area_array = []
    for i in range(0, plot_frame):
        average_area_array.append(average_area[i])
    plt.plot(plot_time_array,average_area_array,'.')
    plt.savefig('average_area_vs_time.png')
    plt.close()

    Lp=[]
    La=[]
    
    Lpf=[]
    Laf=[]
    # the filter here serves the purporse of removing time points that MF has reached the end (if t_prol is non_zero, that is a different story!)
    toggle_mf_filter=False
    MF_init_phase_threshold=50
    for i in range(0,len(dict_from_file['MF position']) ):
        if dict_from_file['MF position'][i]!=0.0:
            Lp+=[dict_from_file['L'][i][0]-dict_from_file['MF position'][i] ]
            La+=[dict_from_file['MF position'][i]-dict_from_file['L'][i][1] ]
            Lpf+=[dict_from_file['L'][i][0]-dict_from_file['MF position'][i] ]
            Laf+=[dict_from_file['MF position'][i]-dict_from_file['L'][i][1] ]
        elif toggle_mf_filter==False:
            if i<MF_init_phase_threshold:
                Lp+=[0.0]
                La+=[dict_from_file['L'][i][0]-dict_from_file['L'][i][1] ]
            else :
                Lp+=[dict_from_file['L'][i][0]-dict_from_file['L'][i][1] ]
                La+=[0.0]
    plt.figure()
    
    time_arrayf=[]
    time_arrayf=time_array[0:len(Lpf)]
    Lp=np.array(Lp)/conversion_r
    La=np.array(La)/conversion_r
    Lpf=np.array(Lpf)/conversion_r
    Laf=np.array(Laf)/conversion_r
    
    plt.ylabel("Posterior Length (μm)")
    plt.xlabel("Time (h)")   
    plt.plot(time_arrayf,Lpf,'.')
    #plt.show()
    plt.savefig('Lp_vs_time.png')
    plt.close()
    
    plt.figure()
    plt.ylabel("Anterior Length (μm)")
    plt.xlabel("Time (h)")
    plt.plot(time_arrayf,Laf,'.')
    #plt.show()                                                                                                                                                                          
    plt.savefig('La_vs_time.png')
    plt.close()

    
    
    
    #Lp=Lp[10:]
    
    coef,rsquared=scatter_and_fit(time_arrayf,Lpf)
    print(coef)
    print("Rsq " + str(rsquared) )
    fitted_MF_speed=coef[0]
    # units: micrometer/tmech = micrometer/(t_mech * hr/t_mech)
    # converting to microm/hr
    #fitted_MF_speed=fitted_MF_speed/(conversion_t_hr)

    
    out_dict={}
    out_dict['MF speed']=fitted_MF_speed
    out_dict['MF linearity']=rsquared
    
    #plot_area_keys=["tissue area","Posterior area","Anterior area"]
    area=np.array(dict_from_file["tissue area"])/((conversion_r)**2)
    anterior_area=np.array(dict_from_file["Anterior area"])/((conversion_r)**2)
    posterior_area=np.array(dict_from_file["Posterior area"])/((conversion_r)**2)

    out_dict['final cell number']=dict_from_file["cell number"][-1]
    out_dict['final average area']=average_area[-1]
    out_dict['final area']  =area[-1]
    out_dict['area']        =area[0:plot_frame]
    out_dict['posterior_area']=posterior_area[0:plot_frame]
    out_dict['anterior_area']=anterior_area[0:plot_frame]
    
    posterior_cell_area=[]
    anterior_cell_area=[]
    posterior_area_average_sum=0.0
    anterior_area_average_sum=0.0
    non_zero_frame = 0
    for i in range(0,len(posterior_area)):
        pn=dict_from_file["Posterior cell number"][i]
        if pn==0:
                posterior_cell_area+=[0]
        else:
                posterior_cell_area+=[posterior_area[i]/pn]
                
        an=dict_from_file["Anterior cell number"][i]
        if an==0:
                anterior_cell_area+=[0]
        else:
                anterior_cell_area+=[anterior_area[i]/an]
        if an != 0 and pn != 0:
            posterior_area_average_sum+=posterior_area[i]/pn
            anterior_area_average_sum+=anterior_area[i]/an
            non_zero_frame += 1
            
    out_dict['posterior_cell_area']=posterior_cell_area[0:plot_frame]
    out_dict['anterior_cell_area']=anterior_cell_area[0:plot_frame]
    # This average is not perfect
    # it includes zeroes, where proliferation ended, and simulation static result. 
    # perhaps we need to stop the simulation when MF reaches anterior to avoid this
    average_posterior_cell_area = posterior_area_average_sum/non_zero_frame
    average_anterior_cell_area = anterior_area_average_sum/non_zero_frame
    out_dict['average_posterior_cell_area']=average_posterior_cell_area
    out_dict['average_anterior_cell_area']=average_anterior_cell_area
    
    total_area_jump = []
    posterior_area_jump = []
    anterior_area_jump = []
    average_total_area_jump = []
    average_posterior_area_jump = []
    average_anterior_area_jump = []
    for i in range(0, plot_frame-1):
        total_area_jump.append(out_dict['area'][i+1]-out_dict['area'][i])
        posterior_area_jump.append(out_dict['posterior_area'][i+1]-out_dict['posterior_area'][i])
        anterior_area_jump.append(out_dict['anterior_area'][i+1]-out_dict['anterior_area'][i])
        average_total_area_jump.append(average_area_array[i+1]-average_area_array[i])
        average_posterior_area_jump.append(out_dict['posterior_cell_area'][i+1]-out_dict['posterior_cell_area'][i])
        average_anterior_area_jump.append(out_dict['anterior_cell_area'][i+1]-out_dict['anterior_cell_area'][i])
    
    plot_time_array_jump = plot_time_array[0:plot_frame-1]
    yname="area_jump"
    plt.figure()
    plt.ylabel(yname+" ($μm^2$)")
    plt.xlabel("Time (h)")
    total_area_label, = plt.plot(plot_time_array_jump,total_area_jump,color = 'green',marker='.', label='total area')
    posterior_area_label, =plt.plot(plot_time_array_jump,posterior_area_jump,color = 'blue',marker='.', label='posterior area')
    anterior_area_label, =plt.plot(plot_time_array_jump,anterior_area_jump,color = 'red',marker='.', label='anterior area')
    plt.legend([total_area_label, posterior_area_label, anterior_area_label])
    plt.savefig(yname+"_vs_time.png")
    plt.close()
    
    yname="cell_area_jump"
    plt.figure()
    plt.ylabel(yname+" ($μm^2$)")
    plt.xlabel("Time (h)")
    average_total_area_label, = plt.plot(plot_time_array_jump,average_total_area_jump,color = 'green',marker='.', label='total area')
    average_posterior_area_label, =plt.plot(plot_time_array_jump,average_posterior_area_jump,color = 'blue',marker='.', label='posterior area')
    average_anterior_area_label, =plt.plot(plot_time_array_jump,average_anterior_area_jump,color = 'red',marker='.', label='anterior area')
    plt.legend([average_total_area_label, average_posterior_area_label, average_anterior_area_label])
    plt.savefig(yname+"_vs_time.png")
    plt.close()

    
    ave_num_sides_MF = []
    for i in range (0, len(dict_from_file["average_number_of_sides_in_MF"])):
        if dict_from_file["average_number_of_sides_in_MF"][i] != 0.0 and np.isnan(dict_from_file["average_number_of_sides_in_MF"][i]) != True:
            ave_num_sides_MF += [dict_from_file["average_number_of_sides_in_MF"][i]]
    out_dict['average_number_of_sides_in_MF']=ave_num_sides_MF

    ave_area_MF = []

    for i in range(0, len(dict_from_file["average_area_in_MF"])):
        if dict_from_file["average_area_in_MF"][i] != 0.0 and np.isnan(dict_from_file["average_area_in_MF"][i]) != True:
            ave_area_MF += [ dict_from_file["average_area_in_MF"][i]/((conversion_r)**2) ]

    out_dict['average_area_in_MF'] = ave_area_MF

    
    MF_shape = []
    for i in range(0, len(dict_from_file["MF_shape"])):
        if dict_from_file["MF_shape"][i] != 0.0 and np.isnan(dict_from_file["MF_shape"][i]) != True:
                MF_shape += [ dict_from_file["MF_shape"][i] ]
    out_dict['MF_shape'] = MF_shape

    out_dict['parameters']=parameters
    MF_last=dict_from_file['MF position'][-1]
    
    out_dict['Lp']=Lpf
    out_dict['La']=Laf
    
    out_dict["Posterior cell number"]=dict_from_file["Posterior cell number"][0:plot_frame]
    out_dict["Anterior cell number"]=dict_from_file["Anterior cell number"][0:plot_frame]
    out_dict["cell number"]=dict_from_file["cell number"][0:plot_frame]
    out_dict["cell_death"]=dict_from_file["cell_death"][0:plot_frame]
    out_dict["cell_division"]=dict_from_file["cell_division"][0:plot_frame]
    
    out_dict["mitotic_frequency"] = mitotic_frequency_list
    out_dict["normalized_time_list"] = normalized_time_list
    out_dict["mitotic_number_list"] = mitotic_number_list
    out_dict["total_mitotic_number"] = total_mitotic_number
    
    if MF_last==0.0:
        out_dict['finished']=True
    else:
        out_dict['finished']=False
    return out_dict




#print(analyze())
