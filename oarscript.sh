#!/bin/bash            
#OAR -l /nodes=8/core=8, walltime=8:00:00
#OAR -t frag_ok
#OAR --project icechrono

echo OAR initialized
#!/usr/bin/env nix-shell                                                    
#!nix-shell -I /scratch/nix /home/jcb232/nix-tests/python/shell.drv -i bash 
                                                            
echo nix-shell activated
               
source /home/jcb232/nix-tests/python/venv/bin/activate      

echo virtual environment activated

cp -r /home/jcb232/IceChrono_new/AICC2012-VLR /scratch/PROJECTS/pr-icechrono/

echo i-o directory copied, running IceChrono

which python
which mpirun

mpirun --machinefile $OAR_NODE_FILE -mca plm_rsh_agent "oarsh" python $OAR_WORKING_DIRECTORY/IceChrono.py /scratch/PROJECTS/pr-icechrono/AICC2012-VLR/

echo IceChrono terminated, copying directory

cp -r /scratch/PROJECTS/pr-icechrono/AICC2012-VLR /home/jcb232/IceChrono_new/

echo Directory copied

rm -r /scratch/PROJECTS/pr-icechrono/AICC2012-VLR

echo Scratch directory deleted, job done.
