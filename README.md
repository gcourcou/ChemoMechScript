# ChemoMechScript

This repository contains the necessarily materials to reproduce the results in "Avalanches During Epithelial Tissue Growth; Uniform Growth and a Drosophila Eye Disc Model" by Courcoubetis et. al. 


There are two models herein: A uniform growth model ran by executing script_uniform_growth.py and an Drosophila eye disc model by executing script_bcsource.py.

For explanation of parameters, see parameter explanation.txt file. 

--------------------------------------------------------------------------------------------------------------------------------- 

To run this code please install the tyssue package dependencies: conda install --only-deps -c conda-forge tyssue

Then, clone the github repository  https://github.com/gcourcou/tyssue.

While in the master branch, install the tyssue package in the cloned directory: python setup.py install

Note the path of the installation and add it to the sys.argv path in the scripts to run our custom tyssue. 

Note : In future releases of tyssue, our code will likely be merged so that one can simply install tyssue and proceed.

After completing the installation, clone the current repository and execute either of the two scripts to run the simulation.

E.g. python script_uniform_growth.py example_name will execute eye disc growth model with the output saved on a directory named example_name .

The remaining files (e.g. Plotter_mate.py, single_dir_analysis.py, cross_dir_analysis.py, import_analyze_script_out.py, mitotic_frequency.py ) are used to analyze the output of the simulations and produce plots of the results.

Note: Some files such as run_george/run_chi, master.py and gen2.py were used to submit jobs using slurm on a remote machine and can be ignored by third parties.
