#job name: {NAME}
#this script's location: {SCRIPT_LOCATION}
#!/bin/bash
#PBS -N {NAME}
#PBS -l nodes=1:ppn=14
#PBS -l walltime=5:00
#PBS -q pace-ice
#PBS -o output_dft

source /nv/pace-ice/bcomer3/ug_esp_env.sh

{COMMAND}