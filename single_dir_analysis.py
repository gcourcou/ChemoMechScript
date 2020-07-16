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

#print(targets)
dir_in_question=str(input('please enter name of target dir' ))
#os.chdir(top_dir+directory)
temp=analyze(bottom=dir_in_question)
print(temp)

