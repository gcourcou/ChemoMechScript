import os
os.system("for d in */; do cp run_george_gen2 \"$d\"; done")
os.system("for d in ./*/ ; do (cd \"$d\" && sbatch run_george_gen2); done")
