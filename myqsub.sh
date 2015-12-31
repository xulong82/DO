#!/bin/sh

FILE=$1

qsub -l nodes=1:ppn=32,walltime=59:59:59 $FILE 

