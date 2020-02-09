#!/bin/sh
export HARLEM_HOME="/share/apps/HARLEM"
export PYTHONPATH="/share/apps/InterX_Simulations/Python3.7:$PYTHONPATH"
export LD_LIBRARY_PATH="/share/apps/anaconda3/lib/python3.7/site-packages/wx:$LD_LIBRARY_PATH"
#export LD_LIBRARY_PATH="/share/apps/intel/mkl/lib/intel64_LIN:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH"
python3 -m harlempy.start_harlem 
