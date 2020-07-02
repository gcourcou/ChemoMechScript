import shutil
import os
import sys
dir_N=int(input('please enter number of dirs' ))

parameter_range=[3*i for i in range(0,dir_N)]
for item in parameter_range :
    try : 
        os.mkdir(str(item)) 
    except:
        keep_old=bool(int(input('keep_old? 0False 1True dirN ' + str(item) )))
        if keep_old==True:
            print("exists kept old")
        elif keep_old==False:
            shutil.rmtree(str(item))
            os.mkdir(str(item))
        else:
            print("input eror")

os.system("for d in */; do cp parameters.txt \"$d\"; done")
os.system("for d in */; do cp run_george \"$d\"; done")
os.system("for d in */; do cp script_bcsource.py \"$d\"; done")
os.system("for d in */; do cp division_functions_aegerter.py \"$d\"; done")
os.system("for d in */; do cp realistic_tissue.hf5 \"$d\"; done")
# copy a file in al lsub dirs for d in */; do cp water.txt "$d"; done
# for d in ./*/ ; do (cd "$d" && somecommand); done do a command in all subdirs

mod_str=str(input('please enter name of target parameter (eg y_dif)' ))
for root, dirs, files in os.walk(".", topdown=True):
    print(dirs)
    for name in dirs:
        os.chdir(os.path.join(root, name))
        cwd=os.getcwd()
        cwd=cwd.split('/')[-1]
        print(cwd)
        magnitude=float(cwd)
        File_object=open("parameters.txt", "r")                                               
        out=open("mod_parameters.txt", "w")
        for aline in File_object:
            values = aline.split()
            #print(values)
            s=" "
            #mod_str=str(input('please enter name of target parameter (eg y_dif)' ))
            if values[0]==mod_str:
                values[1]=str(float(values[1])*(magnitude/10.))
                theline=s.join(values)
                out.write(theline+"  \n")
                #print(theline)
            else:
                out.write(aline)
        out.close()
        File_object.close()
        os.chdir("../")

os.system("for d in ./*/ ; do (cd \"$d\" && mv mod_parameters.txt parameters.txt); done")
os.system("for d in ./*/ ; do (cd \"$d\" && sbatch run_george); done")
