import matplotlib.pyplot as plt
import numpy as np
np.pi

def  master_P_A_area(file="script_out.txt",division_type="uniform") :
      with open(file) as f:  
         dict_from_file = eval(f.read()) 



     plt.figure(figsize=(10,5),dpi=400)

     plt.plot(dict_from_file["Posterior area"])

     plt.title("Posterior area")
     plt.xlabel("time",fontsize=20)
     plt.ylabel("P_area",fontsize=20)
     plt.savefig("P_area"+str(division_type)+".png",bbox_inches='tight',dpi=400)
     plt.close()

     plt.cla()

     plt.plot(dict_from_file["Anterior area"])

     plt.title("Anterior area")
     plt.xlabel("time",fontsize=20)
     plt.ylabel("A_area",fontsize=20)
     plt.savefig("A_area"+str(division_type)+".png",bbox_inches='tight',dpi=400)
     plt.close()
