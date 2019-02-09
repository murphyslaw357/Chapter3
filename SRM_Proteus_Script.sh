#!/bin/bash
#$ -S /bin/bash
#$ -M srm329@drexel.edu
#$ -m beas
#$ -j y
#$ -cwd
#$ -P nieburPrj
#$ -pe fixed16 16
#$ -l h_rt=12:00:00
#$ -l m_mem_free=7G
#$ -l h_vmem=63G
#$ -l matlab=1

. /etc/profile.d/modules.sh
module load shared
module load sge/univa
module load proteus
module load matlab
cd '/mnt/HA/groups/nieburGrp/Shaun/NewtonRaphsonHeatBalance/'

matlab -nodisplay -nodesktop -nosplash -noFigureWindows -r "run ConvergenceProof.m"