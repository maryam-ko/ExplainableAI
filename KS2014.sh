#!/bin/bash
#$ -cwd              # Set the working directory for the job to the current directory
#$ -pe smp 4         # Request 8 cores (adjust based on your needs)
#$ -l h_rt=2:0:0     # Request 2 hours of runtime (adjust as needed)
#$ -l h_vmem=12G     # Request 12 GB of RAM
#$ -o KS2014.csv
#$ -e KS2014_error.txt

# Load the Python module
module load python/3.11.7-gcc-12.2.0

# Activate your virtual environment
source ~/myenv/bin/activate

python /data/home/bt24990/maryam-ko-QMUL-MSc-Project/01_input_data/preprocessing_scripts/KS2014.py > KS2014.csv
