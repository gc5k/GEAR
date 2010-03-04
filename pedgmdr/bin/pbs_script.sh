#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o Test_PBS
#PBS -j oe
#PBS -m e 
#PBS -M gc5k@virginia.edu

cd $PBS_O_WORKDIR
module java/1.5
java -jar null.jar

