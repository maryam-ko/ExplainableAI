#!/bin/bash
#$ -cwd              
#$ -pe smp 1         
#$ -l h_rt=24:0:0     
#$ -l h_vmem=2G 
#$ -o /data/home/bt24990/ExplainableAI/job_output
#$ -j y
#$ -N LRsubset150
#$ -M m.koddus@se24.qmul.ac.uk   
#$ -m e

module load python/3.11.7-gcc-12.2.0
source ~/myenv/bin/activate

python /data/home/bt24990/ExplainableAI/06_models/linear_regression/062_linear_regression_subset_coefficients.py