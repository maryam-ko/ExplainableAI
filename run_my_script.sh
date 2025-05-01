#!/bin/bash
#$ -cwd              # Set the working directory for the job to the current directory
#$ -pe smp 4         # Request 8 cores (adjust based on your needs)
#$ -l h_rt=2:0:0     # Request 2 hours of runtime (adjust as needed)
#$ -l h_vmem=12G     # Request 12 GB of RAM
#$ -o HH2022.csv     # Redirect standard output to HH2022.csv
#$ -e HH2022_error.txt  # Redirect error messages to HH2022_error.txt

module load python/3.11.7-gcc-12.2.0
source ~/myenv/bin/activate

python /data/home/bt24990/maryam-ko-QMUL-MSc-Project/01_input_data/preprocessing_scripts/HH2022.py > HH2022.csv
