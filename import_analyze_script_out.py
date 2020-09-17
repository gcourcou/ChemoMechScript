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
    plot_area_keys=["tissue area","Posterior area","Anterior area"]
    #names used for plots
    plot_area_names=["Total area","Posterior area","Anterior area"]
    plot_number_of_sides_in_MF=["average_number_of_sides_in_MF", "MF_shape"]
    plot_area_in_MF=["average_area_in_MF"]
    for i in range(0,len(plot_area_keys)):
        key=plot_area_keys[i]
        yname=plot_area_names[i]
        plt.figure()
        plt.ylabel(yname+" ($μm^2$)")
        plt.xlabel("Time (h)")
        area=np.array(dict_from_file[key])/((conversion_r)**2)
        plt.plot(time_array,area,'.')
        plt.savefig(yname+"_vs_time.png")
        plt.close()
    for i in range(0,len(plot_number_of_sides_in_MF)):
        key=plot_number_of_sides_in_MF[i]
        yname="average_number_of_sides_in_MF"
        plt.figure()
        plt.ylabel(yname+" ($μm^2$)")
        plt.xlabel("Time (h)")
        ylist = []
        time_array_1 = []
        for i in range(0, len(time_array)):
            if dict_from_file[key][i] != 0.0 and np.isnan(dict_from_file[key][i]) != True:
                ylist += [dict_from_file[key][i]]
                time_array_1 += [time_array[i]]
        plt.plot(time_array_1,ylist,'.')
        plt.savefig(yname+"_vs_time.png")
        plt.close()
    for i in range(0,len(plot_area_in_MF)):
        key=plot_area_in_MF[i]
        yname="average_area_in_MF"
        plt.figure()
        plt.ylabel(yname+" ($μm^2$)")
        plt.xlabel("Time (h)")
        area= []
        time_array_1 = []
        for i in range(0, len(time_array)):
            if dict_from_file[key][i] != 0.0 and np.isnan(dict_from_file[key][i]) != True:
                area += [dict_from_file[key][i]/((conversion_r)**2)]
                time_array_1 += [time_array[i]]
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
    plt.plot(time_array,dict_from_file["cell number"],'.')
    #plt.show()
    plt.savefig('cell_number_vs_time.png')
    plt.close()
    
    average_area=np.array(dict_from_file["tissue area"])/np.array(dict_from_file["cell number"])
    average_area=average_area/((conversion_r)**2)
    plt.figure()
    plt.ylabel("Average cell area ($μm^2$)")
    plt.xlabel("Time (h)")
    plt.plot(time_array,average_area,'.')
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
    plt.plot(time_array,Lp,'.')
    #plt.show()
    plt.savefig('Lp_vs_time.png')
    plt.close()
    
    plt.figure()
    plt.ylabel("Anterior Length (μm)")
    plt.xlabel("Time (h)")
    plt.plot(time_array,La,'.')
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


    out_dict['final average area']=average_area[-1]
    out_dict['final area']  =area[-1]
    out_dict['area']        =area
    out_dict['posterior_area']=posterior_area
    out_dict['anterior_area']=anterior_area
    
    posterior_cell_area=[]
    anterior_cell_area=[]
    for i in range(0,len(posterior_area)):
        cn=dict_from_file["Posterior cell number"][i]
        if cn==0:
                posterior_cell_area+=[0]
        else:
                posterior_cell_area+=[posterior_area[i]/cn]
        cn=dict_from_file["Anterior cell number"][i]
        if cn==0:
                anterior_cell_area+=[0]
        else:
                anterior_cell_area+=[anterior_area[i]/cn]
    out_dict['posterior_cell_area']=posterior_cell_area
    out_dict['anterior_cell_area']=anterior_cell_area
    # This average is not perfect
    # it includes zeroes, where proliferation ended, and simulation static result. 
    # perhaps we need to stop the simulation when MF reaches anterior to avoid this
    out_dict['average_posterior_cell_area']=np.average(posterior_cell_area)
    out_dict['average_anterior_cell_area']=np.average(anterior_cell_area)
    
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
    
    out_dict['Lp']=Lp
    out_dict['La']=La
    out_dict["Posterior cell number"]=dict_from_file["Posterior cell number"]
    out_dict["Anterior cell number"]=dict_from_file["Anterior cell number"]
    out_dict["cell number"]=dict_from_file["cell number"]

    if MF_last==0.0:
        out_dict['finished']=True
    else:
        out_dict['finished']=False
    return out_dict




#print(analyze())
