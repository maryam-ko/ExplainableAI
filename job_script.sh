#!/bin/bash
#$ -cwd              
#$ -pe smp 4         
#$ -l h_rt=1:0:0     
#$ -l h_vmem=12G     
#$ -o FS2019.csv     
#$ -e FS2019_error.txt  
#$ -M m.koddus@se24.qmul.ac.uk   
#$ -m e                        

module load python/3.11.7-gcc-12.2.0
source ~/myenv/bin/activate

python /data/home/bt24990/maryam-ko-QMUL-MSc-Project/01_input_data/preprocessing_scripts/FS2019.py
 