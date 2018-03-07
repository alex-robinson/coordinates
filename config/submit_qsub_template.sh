#!/bin/bash
#$ -V                             # Ensure user enivronment variables are available
#$ -cwd                           # To use the current directory
#$ -m ae                          # Send mail when job is aborted (a), begins (b) and ends (e)
#$ -N  yelmo                      # (nombre del trabajo)
#$ -o ./out.out                   # (fichero de salida)   $JOB_ID
#$ -e ./out.err                   # (fichero de error)


### jalv: this environnement variable limits the number of CPUs used in eolo (=1, means one cpu)
###export OMP_NUM_THREADS=1
# Currently in .bashrc, but maybe should be done here.

# Run the job
./libyelmo/bin/yelmo_eismint.x 

