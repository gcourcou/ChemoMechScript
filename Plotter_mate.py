#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:13:19 2020

@author: georgecourcoubetis
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_picture(axis,f_name,f_format,label="A$)$",a_spect="equal"):
    bn=plt.imread(f_name,format=f_format)
    axis.imshow(bn,aspect=a_spect)
    axis.text(-0.02, 0.5, label, transform=axis.transAxes,
              fontsize=20, fontweight='bold', va='top', ha='right')
    axis.axis('off')


# Initialize pictures
s_trimmed=True
#plt.figure(figsize=(20,20))

n_rows=2
n_columnns=5


fig, ax = plt.subplots(nrows=n_rows, ncols=n_columnns, figsize=(15,3*n_rows),dpi=400)
##,gridspec_kw={'width_ratios': [1.4,1,1]}


#1
#plot_picture(ax[0][0],'spot_number.png','png',label="A$)$")
#plot_picture(ax[0][1],'filtered_pixelval0090-1.tif','tiff',label="B$)$")
#plot_picture(ax[0][2],'density_histogram_2d_FOI.png','png',label="C$)$")
#plot_picture(ax[0][3],'filtered_pixelval0090_zoomed_80_80_100-1.tiff','tiff',label="D$)$")
dirs=[0.1,0.2,0.30000000000000004,0.4,0.5,0.6,0.7000000000000001,0.8,0.9,1.0]

plot_string=str(input("enter plot name eg. imageareavsxarea06.png"))
#plot_string="imageareavsxarea06.png"
dir_string=str(input("dir string eg. /out_0/"))
#dir_string="/out_0/"
c=0
for key in dirs:
    if c<5:
        plot_picture(ax[0][c%5],str(key)+dir_string+plot_string,'png',label=str(key)[0:3])
    else:
        plot_picture(ax[1][c%5],str(key)+dir_string+plot_string,'png',label=str(key)[0:3])
    c+=1
#2
#plot_picture(ax[0][0],'area90.png','png',label="A$)$")




fig.tight_layout()
plt.savefig('test.pdf', bbox_inches='tight')
#plt.show()
#plt.close()
