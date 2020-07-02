import matplotlib
matplotlib.use("Agg")

import matplotlib.pylab as plt     
import numpy as np

from scipy.stats import chisquare

from sklearn.metrics import r2_score

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
def  scatter_and_fit(target_y,text_x="t_mech",text_y="Lp",text="none"):
    #fig=plt.figure(figsize=(7,5),dpi=calc)
    plt.figure()
    ax=plt.axes()
    target_x=np.linspace(1,len(target_y), num=len(target_y))-1
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
    t_mech=1
    
    proliferation_magnitude = float(parameters["proliferation_magnitude"])
    
    with open(bottom+"script_out.txt") as f:  
        dict_from_file = eval(f.read()) 
    
    plt.figure()
    plt.ylabel("Area (mm^2)")
    plt.xlabel("t_mech")
    conversion_r=0.23232956642491454
    #pixel/microm
    conversion_t=915.3565322296259
    #s/t_mech
    conversion_t_hr=conversion_t*0.000277778
    #hr/t_mech
    area=np.array(dict_from_file["tissue area"])/((conversion_r)**2)
    plt.plot(area,'.')
    plt.savefig('area_vs_time.png')
    plt.close()
    Lp=[]
    La=[]
    # the filter here serves the purporse of removing time points that MF has reached the end (if t_prol is non_zero, that is a different story!)
    toggle_mf_filter=True
    for i in range(0,len(dict_from_file['MF position']) ):
        if dict_from_file['MF position'][i]!=0.0:
                Lp+=[dict_from_file['L'][i][0]-dict_from_file['MF position'][i] ]
                La+=[dict_from_file['MF position'][i]-dict_from_file['L'][i][1] ]
        elif toggle_mf_filter==False:
                Lp+=[dict_from_file['L'][i][0]-dict_from_file['MF position'][i] ]
                La+=[dict_from_file['MF position'][i]-dict_from_file['L'][i][1] ]
    plt.figure()
    Lp=np.array(Lp)/conversion_r
    La=np.array(La)/conversion_r
    plt.ylabel("Lp")
    plt.xlabel("t_mech")   
    plt.plot(Lp,'.')
    #plt.show()
    plt.savefig('Lp_vs_time.png')
    plt.close()
    
    plt.figure()
    plt.ylabel("La")
    plt.xlabel("t_mech")
    plt.plot(La,'.')
    #plt.show()                                                                                                                                                                          
    plt.savefig('La_vs_time.png')
    plt.close()

    
    
    
    #Lp=Lp[10:]
    
    coef,rsquared=scatter_and_fit(Lp)
    print(coef)
    print("Rsq " + str(rsquared) )
    fitted_MF_speed=coef[0]
    # units: micrometer/tmech = micrometer/(t_mech * hr/t_mech)
    # converting to microm/hr
    fitted_MF_speed=fitted_MF_speed/(conversion_t_hr)


    average_area=np.array(dict_from_file["tissue area"])/np.array(dict_from_file["cell number"])
    plt.figure()
    plt.ylabel("Average Cell Area")
    plt.xlabel("t_mech")
    plt.plot(average_area,'.')
    plt.savefig('average_area_vs_time.png')
    plt.close()

    
    out_dict={}
    out_dict['MF speed']=fitted_MF_speed
    out_dict['MF linearity']=rsquared
    out_dict['final average area']=average_area[-1]
    out_dict['final area']  =area[-1]
    out_dict['area']        =area
    out_dict['parameters']=parameters
    MF_last=dict_from_file['MF position'][-1]
    
    out_dict['Lp']=Lp
    out_dict['La']=La
    if MF_last==0.0:
        out_dict['finished']=True
    else:
        out_dict['finished']=False
    return out_dict




#print(analyze())




















