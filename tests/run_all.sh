#!/bin/sh                                                               
# run  input datafile from Ishara                                     
# then compare results with *_test                                              
# Shujie, 201912                                                                

for filename in  input/*.dat ;do
        ./tstfit16 $filename
        mv ${filename::-4}.out ${filename::-4}.out_shujie
done


for filename in input/*.dat ;do
        diff ${filename::-4}.out_shujie ${filename::-4}.out_new >>${filename::-4}.diff
done